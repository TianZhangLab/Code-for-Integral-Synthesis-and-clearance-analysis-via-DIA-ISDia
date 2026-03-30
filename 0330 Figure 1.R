#Figure 1

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(colorspace)
library(scales)
library(patchwork)
library(grid)

DIA <- read.delim(
  "D:/personal/UVA/Data analysis/202412 SILAC/DIA combination/20251207_134746_V196_DIA_all_ratios_Report 1.tsv",
  header = TRUE, sep = "\t"
)

ratio_levels_main <- c("1/27","1/9","1/3","1/1","3/1","9/1","27/1")

core_base <- c(
  "1/27" = "#4DADDE",
  "1/9"  = "#33B6A5",
  "1/3"  = "#54BF7A",
  "1/1"  = "#DE6E5C",
  "3/1"  = "#ffd24a",
  "9/1"  = "#b57db6",
  "27/1" = "#9aa3ad"
)

safe_cv <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 2) return(NA_real_)
  m <- mean(x)
  if (!is.finite(m) || m == 0) return(Inf)
  stats::sd(x) / m * 100
}

parse_ratio_value <- function(r) {
  a <- as.numeric(sub("/.*", "", r))
  b <- as.numeric(sub(".*/", "", r))
  a / b
}

clean_DIA <- function(df) {
  ch1_cols <- grep("\\.PG\\.MS2Channel1$", names(df), value = TRUE)
  ch2_cols <- grep("\\.PG\\.MS2Channel2$", names(df), value = TRUE)
  
  ch1_pref <- sub("\\.PG\\.MS2Channel1$", "", ch1_cols)
  ch2_pref <- sub("\\.PG\\.MS2Channel2$", "", ch2_cols)
  
  ch2_cols_aligned <- ch2_cols[match(ch1_pref, ch2_pref)]
  stopifnot(!any(is.na(ch2_cols_aligned)))
  
  for (i in seq_along(ch1_cols)) {
    c1 <- ch1_cols[i]
    c2 <- ch2_cols_aligned[i]
    idx <- is.na(df[[c1]]) | is.na(df[[c2]])
    df[idx, c(c1, c2)] <- NA_real_
  }
  df
}

prep_dia_wide_lower <- function(df) {
  df %>%
    pivot_longer(
      cols = matches("\\.PG\\.MS2Channel[12]$"),
      names_to = "col",
      values_to = "intensity"
    ) %>%
    mutate(
      channel = if_else(str_detect(col, "MS2Channel1$"), "light", "heavy"),
      Ratio = str_extract(col, "HL_\\d+_\\d+") %>% str_replace("HL_(\\d+)_(\\d+)", "\\1/\\2"),
      sample_key = str_remove(col, "\\.PG\\.MS2Channel[12]$"),
      sample_key = str_remove(sample_key, "^X\\.\\d+\\.\\."),
      intensity = suppressWarnings(as.numeric(intensity)),
      intensity = if_else(is.na(intensity) | intensity <= 0, NA_real_, intensity)
    ) %>%
    group_by(Ratio) %>%
    mutate(rep = match(sample_key, sort(unique(sample_key)))) %>%
    ungroup() %>%
    transmute(Protein.ID = PG.ProteinAccessions, Ratio, rep, channel, intensity) %>%
    pivot_wider(names_from = channel, values_from = intensity)
}

prep_dia_wide_upper <- function(df) {
  df %>%
    pivot_longer(
      cols = matches("\\.PG\\.MS2Channel[12]$"),
      names_to = "col",
      values_to = "Intensity"
    ) %>%
    mutate(
      Channel = if_else(str_detect(col, "MS2Channel1$"), "Light", "Heavy"),
      Ratio = str_extract(col, "HL_\\d+_\\d+") %>% str_replace("HL_(\\d+)_(\\d+)", "\\1/\\2"),
      sample_key = str_remove(col, "\\.PG\\.MS2Channel[12]$"),
      sample_key = str_remove(sample_key, "^X\\.\\d+\\.\\."),
      Intensity = suppressWarnings(as.numeric(Intensity)),
      Intensity = if_else(is.na(Intensity) | Intensity <= 0, NA_real_, Intensity)
    ) %>%
    group_by(Ratio) %>%
    mutate(Rep = match(sample_key, sort(unique(sample_key)))) %>%
    ungroup() %>%
    transmute(Protein.ID = PG.ProteinAccessions, Ratio, Rep, Channel, Intensity) %>%
    pivot_wider(names_from = Channel, values_from = Intensity)
}

prep_dia_long <- function(df) {
  df %>%
    pivot_longer(
      cols = matches("\\.PG\\.MS2Channel[12]$"),
      names_to = "col",
      values_to = "Intensity"
    ) %>%
    mutate(
      Intensity_Type = if_else(str_detect(col, "MS2Channel1$"), "Light", "Heavy"),
      Ratio = str_extract(col, "HL_\\d+_\\d+") %>% str_replace("HL_(\\d+)_(\\d+)", "\\1/\\2"),
      sample_key = str_remove(col, "\\.PG\\.MS2Channel[12]$"),
      sample_key = str_remove(sample_key, "^X\\.\\d+\\.\\."),
      Intensity = suppressWarnings(as.numeric(Intensity)),
      Intensity = if_else(is.na(Intensity) | Intensity <= 0, NA_real_, Intensity)
    ) %>%
    group_by(Ratio) %>%
    mutate(Rep = match(sample_key, sort(unique(sample_key)))) %>%
    ungroup() %>%
    transmute(Protein.ID = PG.ProteinAccessions, Ratio, Rep, Intensity_Type, Intensity)
}

summarise_ratio_counts <- function(wide_df) {
  wide_df %>%
    mutate(valid = !is.na(light) & !is.na(heavy)) %>%
    group_by(Protein.ID, Ratio) %>%
    summarise(
      n_rep    = sum(valid),
      light_cv = safe_cv(light[valid]),
      heavy_cv = safe_cv(heavy[valid]),
      total_cv = safe_cv((light + heavy)[valid]),
      hl_cv    = safe_cv((heavy / light)[valid]),
      cv_ok = (n_rep >= 3) &
        !is.na(light_cv) & light_cv < 20 &
        !is.na(heavy_cv) & heavy_cv < 20 &
        !is.na(total_cv) & total_cv < 20 &
        !is.na(hl_cv)    & hl_cv < 20,
      .groups = "drop"
    ) %>%
    group_by(Ratio) %>%
    summarise(
      n1   = sum(n_rep >= 1),
      n2   = sum(n_rep >= 2),
      n3   = sum(n_rep >= 3),
      n3cv = sum(cv_ok),
      .groups = "drop"
    )
}

make_plot_data <- function(sum_df, method_name, ratio_levels) {
  sum_df2 <- sum_df %>%
    mutate(Ratio = factor(Ratio, levels = ratio_levels)) %>%
    tidyr::complete(Ratio, fill = list(n1 = 0L, n2 = 0L, n3 = 0L, n3cv = 0L)) %>%
    mutate(
      seg1 = pmax(n1 - n2, 0),
      seg2 = pmax(n2 - n3, 0),
      seg3 = pmax(n3 - n3cv, 0),
      seg4 = pmax(n3cv, 0)
    )
  
  counts <- sum_df2 %>%
    dplyr::select(Ratio, seg1, seg2, seg3, seg4) %>%
    pivot_longer(cols = starts_with("seg"), names_to = "seg", values_to = "Count") %>%
    mutate(
      Occurrence = recode(seg, seg1 = "1", seg2 = "2", seg3 = "3", seg4 = "3(CV<20)"),
      Method = method_name
    ) %>%
    dplyr::select(Ratio, Occurrence, Count, Method)
  
  list(counts = counts)
}

build_fill_values <- function(ratio_levels, core_base) {
  ratio_base <- core_base[ratio_levels]
  
  shade_one_ratio <- function(col) {
    c(
      "1"        = colorspace::lighten(col, 0.55),
      "2"        = colorspace::lighten(col, 0.35),
      "3"        = colorspace::lighten(col, 0.2),
      "3(CV<20)" = col
    )
  }
  
  fill_values <- unlist(lapply(names(ratio_base), function(r) {
    sc <- shade_one_ratio(ratio_base[[r]])
    setNames(sc, paste0(r, "__", names(sc)))
  }))
  
  list(
    fill_values = fill_values,
    occ_levels_legend = c("1", "2", "3", "3(CV<20)"),
    legend_breaks = paste("1/1", c("1","2","3","3(CV<20)"), sep="__"),
    legend_greys = c("grey85","grey65","grey55","grey35")
  )
}

plot_occurrence_bar_dia <- function(DIA_clean, ratio_levels, core_base) {
  dia_wide <- prep_dia_wide_lower(DIA_clean)
  dia_sum  <- summarise_ratio_counts(dia_wide)
  dia_pd   <- make_plot_data(dia_sum, "DIA", ratio_levels)
  
  plot_df <- dia_pd$counts %>%
    mutate(
      Ratio = factor(Ratio, levels = ratio_levels),
      Occurrence = factor(Occurrence, levels = c("3(CV<20)","3","2","1")),
      fill_id = paste(Ratio, Occurrence, sep = "__")
    )
  
  seg_label_df <- plot_df %>%
    arrange(Ratio, Occurrence) %>%
    group_by(Ratio) %>%
    mutate(
      cum = cumsum(Count),
      y_mid = cum - Count / 2,
      seg_label = if_else(Count > 0, as.character(cum), "")
    ) %>%
    ungroup()
  
  fills <- build_fill_values(ratio_levels, core_base)
  
  ggplot(plot_df, aes(x = Ratio, y = Count, fill = fill_id)) +
    geom_col(width = 0.8) +
    geom_text(
      data = seg_label_df,
      aes(x = Ratio, y = y_mid, label = seg_label),
      inherit.aes = FALSE,
      size = 3
    ) +
    scale_fill_manual(
      values = fills$fill_values,
      breaks = fills$legend_breaks,
      labels = fills$occ_levels_legend,
      name = "Occurrence",
      guide = guide_legend(override.aes = list(fill = fills$legend_greys))
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(
      breaks = function(x) sort(unique(c(scales::pretty_breaks()(x), 10000))),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = "H/L Ratio", y = "Protein count") +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = NA, color = NA),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.18, "cm"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    )
}

cv_long <- function(wide_df, method_name = "DIA") {
  wide_df %>%
    mutate(
      Total = Light + Heavy,
      HL = Heavy / Light,
      valid = is.finite(Light) & is.finite(Heavy) & !is.na(Light) & !is.na(Heavy)
    ) %>%
    filter(valid) %>%
    group_by(Protein.ID, Ratio) %>%
    summarise(
      CV_Light = safe_cv(Light),
      CV_Heavy = safe_cv(Heavy),
      Total_intensity_CV = safe_cv(Total),
      HtoL_ratio_CV = safe_cv(HL),
      .groups = "drop"
    ) %>%
    mutate(Acquisition = method_name) %>%
    pivot_longer(
      cols = c(CV_Light, CV_Heavy, Total_intensity_CV, HtoL_ratio_CV),
      names_to = "Group",
      values_to = "CV"
    )
}

make_one_panel <- function(df, group_key, title_name, ratio_levels_main, dia_color = "#9ccea7") {
  data_subset <- df %>%
    filter(Group == group_key) %>%
    mutate(Ratio = factor(Ratio, levels = ratio_levels_main))
  
  ggplot(data_subset, aes(x = Ratio, y = CV)) +
    geom_boxplot(width = 0.75, color = "black", fill = dia_color, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 18, size = 2.6) +
    geom_hline(yintercept = 20, linetype = "dashed", linewidth = 0.6) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_x_discrete(drop = FALSE) +
    labs(title = title_name, x = "H/L ratio", y = "CV (%)") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.18, "cm"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = "none",
      strip.background = element_blank()
    )
}

make_dia_summary_from_wide <- function(dia_wide) {
  dia_wide %>%
    mutate(
      valid = is.finite(Light) & is.finite(Heavy) & !is.na(Light) & !is.na(Heavy) & Light > 0 & Heavy > 0,
      HL = Heavy / Light,
      Total = Heavy + Light
    ) %>%
    filter(valid) %>%
    group_by(Protein.ID, Ratio) %>%
    summarise(
      Total_mean = mean(Total, na.rm = TRUE),
      HtoL_ratio_mean = mean(HL, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Acquisition = "DIA")
}


plot_scatter_with_density_dia_total <- function(df_summary,
                                                ratio_levels = NULL,
                                                color_palette = NULL,
                                                x_lim = c(5, 23),
                                                y_lim = c(-10, 10),
                                                widths = c(1.6, 0.6)) {
  if (is.null(ratio_levels)) ratio_levels <- sort(unique(df_summary$Ratio))
  if (is.null(color_palette)) color_palette <- setNames(rep("#ABD9EF", length(ratio_levels)), ratio_levels)
  
  plot_df <- df_summary %>%
    mutate(
      Ratio = factor(Ratio, levels = ratio_levels),
      log2_HL = log2(HtoL_ratio_mean),
      log2_Total = log2(Total_mean)
    ) %>%
    filter(is.finite(log2_HL), is.finite(log2_Total)) %>%
    filter(log2_HL >= y_lim[1] & log2_HL <= y_lim[2]) %>%
    filter(log2_Total >= x_lim[1] & log2_Total <= x_lim[2])
  
  ratio_line_df <- data.frame(
    Ratio = factor(ratio_levels, levels = ratio_levels),
    y = log2(parse_ratio_value(ratio_levels))
  )
  
  # stats for labels
  sd_data <- plot_df %>%
    group_by(Ratio) %>%
    summarise(
      sd_value = sd(log2_HL, na.rm = TRUE),
      mean_log2_HL = mean(log2_HL, na.rm = TRUE),
      median_log2_HL = median(log2_HL, na.rm = TRUE),
      .groups = "drop"
    )
  
  # positions for median labels on scatter (near right edge, inside panel)
  x_range <- diff(x_lim)
  x_median_scatter <- x_lim[2] - 0.02 * x_range      # number column
  x_median_header  <- x_lim[2] - 0.12 * x_range      # header a bit left to avoid clipping
  y_header_scatter <- min(y_lim[2] - 0.4, max(sd_data$median_log2_HL, na.rm = TRUE) + 0.8)
  
  main_plot <- ggplot(plot_df, aes(x = log2_Total, y = log2_HL, color = Ratio)) +
    geom_point(alpha = 0.5, size = 0.2) +
    coord_cartesian(xlim = x_lim, ylim = y_lim, clip = "off") +
    scale_color_manual(values = color_palette) +
    geom_hline(
      data = ratio_line_df,
      aes(yintercept = y),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    # median values shown at the right edge of scatter panel
    geom_text(
      data = sd_data,
      aes(x = x_median_scatter, y = median_log2_HL, label = round(median_log2_HL, 2)),
      inherit.aes = FALSE,
      color = "black",
      size = 3,
      hjust = 1
    ) +
    annotate(
      "text",
      x = x_median_header,
      y = y_header_scatter,
      label = "median",
      hjust = 1.2,
      size = 4,
      color = "black"
    ) +
    labs(x = "log2(H+L mean intensity)", y = "log2(H/L ratio)", color = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "none",
      # give a bit of room in case text is close to the panel edge
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
  
  # density x-range for s.d. labels on the right side of density panel
  dens_vec <- plot_df$log2_HL[is.finite(plot_df$log2_HL)]
  dens_obj <- if (length(dens_vec) < 2) list(y = 1) else density(dens_vec)
  
  max_x <- max(dens_obj$y)
  x_pos <- max_x * 6.5
  x_text <- x_pos * 1.5
  x_max <- max(1, x_text * 1.2)
  
  density_plot <- ggplot(plot_df, aes(y = log2_HL, color = Ratio)) +
    geom_density(linewidth = 1, fill = NA) +
    scale_color_manual(values = color_palette) +
    # SD labels back on density panel right side
    geom_text(
      data = sd_data,
      aes(x = x_text, y = mean_log2_HL, label = round(sd_value, 2)),
      inherit.aes = FALSE,
      color = "black",
      size = 3,
      hjust = 0
    ) +
    annotate(
      "text",
      x = x_text,
      y = max(sd_data$mean_log2_HL, na.rm = TRUE) + 1,
      label = "s.d.",
      hjust = 0,
      size = 4,
      color = "black"
    ) +
    geom_hline(
      data = ratio_line_df,
      aes(yintercept = y),
      linetype = "dashed",
      color = "gray",
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = c(0, 0.5, 1, 1.5),
      minor_breaks = NULL,
      limits = c(0, x_max)
    ) +
    labs(x = "Density", y = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = y_lim)
  
  main_plot + density_plot +
    plot_layout(widths = widths) +
    plot_annotation(
      title = "DIA",
      theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
    )
}


filter_cv20_keep3 <- function(long_df) {
  wide <- long_df %>%
    pivot_wider(names_from = Intensity_Type, values_from = Intensity)
  
  metrics <- wide %>%
    mutate(
      Total = Light + Heavy,
      HL = Heavy / Light,
      valid = is.finite(Light) & is.finite(Heavy) & !is.na(Light) & !is.na(Heavy)
    ) %>%
    group_by(Protein.ID, Ratio) %>%
    summarise(
      n_rep = sum(valid),
      cv_light = safe_cv(Light[valid]),
      cv_heavy = safe_cv(Heavy[valid]),
      cv_total = safe_cv(Total[valid]),
      cv_hl = safe_cv(HL[valid]),
      pass = (n_rep >= 3) &
        !is.na(cv_light) & cv_light < 20 &
        !is.na(cv_heavy) & cv_heavy < 20 &
        !is.na(cv_total) & cv_total < 20 &
        !is.na(cv_hl) & cv_hl < 20,
      .groups = "drop"
    )
  
  keep_pairs <- metrics %>% filter(pass) %>% dplyr::select(Protein.ID, Ratio)
  
  long_df %>% inner_join(keep_pairs, by = c("Protein.ID", "Ratio"))
}


filter_nrep3 <- function(long_df) {
  wide <- long_df %>%
    pivot_wider(names_from = Intensity_Type, values_from = Intensity)
  
  metrics <- wide %>%
    mutate(
      valid = is.finite(Light) & is.finite(Heavy) & !is.na(Light) & !is.na(Heavy)
    ) %>%
    group_by(Protein.ID, Ratio) %>%
    summarise(
      n_rep = sum(valid),
      pass = n_rep >= 3,
      .groups = "drop"
    )
  
  keep_pairs <- metrics %>%
    filter(pass) %>%
    dplyr::select(Protein.ID, Ratio)
  
  long_df %>%
    inner_join(keep_pairs, by = c("Protein.ID", "Ratio"))
}


normalize_to_ratio11 <- function(long_df, ratio_ref = "1/1") {
  ref <- long_df %>%
    filter(Ratio == ratio_ref) %>%
    group_by(Protein.ID, Intensity_Type) %>%
    summarise(ref_mean = mean(Intensity, na.rm = TRUE), .groups = "drop")
  
  long_df %>%
    left_join(ref, by = c("Protein.ID", "Intensity_Type")) %>%
    mutate(
      Normalized = Intensity / ref_mean,
      log2Rel = log2(Normalized)
    ) %>%
    filter(is.finite(log2Rel))
}

make_boxplot <- function(df_norm, title_text) {
  df_norm <- df_norm %>%
    mutate(
      Ratio = factor(Ratio, levels = ratio_levels_main),
      Intensity_Type = factor(Intensity_Type, levels = c("Heavy", "Light"))
    )
  
  ggplot(df_norm, aes(x = Ratio, y = log2Rel, fill = Intensity_Type)) +
    geom_hline(yintercept = 0, linewidth = 0.6) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9), width = 0.65) +
    #stat_summary(fun = mean, geom = "point", shape = 18, size = 3,position = position_dodge(width = 0.9), color = "black") +
    stat_summary(
      fun = mean, geom = "text",
      aes(label = format(round(after_stat(y), 2), nsmall = 2)),
      position = position_dodge(width = 0.9),
      vjust = -0.6, size = 3.2, color = "black"
    ) +
    scale_fill_manual(values = c("Heavy" = "#F14951", "Light" = "#66B2FF")) +
    coord_cartesian(ylim = c(-5.5, 2.5)) +
    labs(title = title_text, x = "H/L ratio", y = "log2 relative abundance", fill = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.18, "cm"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
}

DIA_clean <- clean_DIA(DIA)

p_occ_dia <- plot_occurrence_bar_dia(DIA_clean, ratio_levels_main, core_base)

dia_wide_cv <- prep_dia_wide_upper(DIA_clean)
cv_all <- cv_long(dia_wide_cv, "DIA") %>%
  mutate(
    Ratio = factor(Ratio, levels = ratio_levels_main),
    Acquisition = factor(Acquisition, levels = "DIA")
  )

p_light <- make_one_panel(cv_all, "CV_Light", "Light intensity", ratio_levels_main)
p_heavy <- make_one_panel(cv_all, "CV_Heavy", "Heavy intensity", ratio_levels_main)
p_total <- make_one_panel(cv_all, "Total_intensity_CV", "Total intensity", ratio_levels_main)
p_hl    <- make_one_panel(cv_all, "HtoL_ratio_CV", "H/L ratio", ratio_levels_main)

p_cv <- (p_light | p_heavy | p_total | p_hl)

dia_summary <- make_dia_summary_from_wide(dia_wide_cv)
p_scatter <- plot_scatter_with_density_dia_total(
  dia_summary,
  ratio_levels = ratio_levels_main,
  color_palette = core_base[ratio_levels_main],
  x_lim = c(5, 23),
  y_lim = c(-10, 10),
  widths = c(1.6, 0.6)
)

dia_long <- prep_dia_long(DIA_clean)
dia_filt <- filter_nrep3(dia_long)
dia_norm <- normalize_to_ratio11(dia_filt, ratio_ref = "1/1")
p_norm <- make_boxplot(dia_norm, "DIA")

p_occ_dia
p_cv
p_scatter
p_norm
