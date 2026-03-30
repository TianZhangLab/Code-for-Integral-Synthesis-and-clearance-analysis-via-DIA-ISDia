
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)
library(minpack.lm)
library(purrr)
library(readxl)
library(tidyverse)
library(patchwork)


df1 <- read.delim("D:/personal/UVA/Data analysis/202412 SILAC/A012_treated/20250204_132507_A012_SILAC_ALL2_Report with gene.tsv", header = TRUE, sep = "\t")



add_total_and_hl_when_paired <- function(df1) {
  
  light_cols <- names(df1)[str_detect(names(df1), "PG\\.MS2Channel1$")]
  
  for (lc in light_cols) {
    base <- str_remove(lc, "\\.PG\\.MS2Channel1$")
    hc   <- paste0(base, ".PG.MS2Channel2")
    
    if (hc %in% names(df1)) {
      
      total_col <- paste0(base, ".PG.Total")
      hl_col    <- paste0(base, ".PG.HL_ratio")
      
      light <- suppressWarnings(as.numeric(df1[[lc]]))
      heavy <- suppressWarnings(as.numeric(df1[[hc]]))
      
      paired <- !is.na(light) & !is.na(heavy)
      
      df1[[total_col]] <- ifelse(paired, light + heavy, NA_real_)
      
      df1[[hl_col]] <- ifelse(
        paired & light != 0,
        heavy / light,
        NA_real_
      )
    }
  }
  
  df1
}

df1_with_metrics <- add_total_and_hl_when_paired(df1)


calculate_mean_cv_4metrics <- function(df, min_rep = 3) {
  
  id_cols <- c("PG.ProteinAccessions", "PG.Genes")
  
  long_df <- df %>%
    pivot_longer(
      cols = -all_of(id_cols),
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    mutate(
      Intensity = suppressWarnings(as.numeric(Intensity)),
      
      Metric = case_when(
        str_detect(Sample, "PG\\.MS2Channel1$") ~ "Light",
        str_detect(Sample, "PG\\.MS2Channel2$") ~ "Heavy",
        str_detect(Sample, "PG\\.Total$") ~ "Total",
        str_detect(Sample, "PG\\.HL_ratio$") ~ "HL",
        TRUE ~ NA_character_
      ),
      
      Treatment = case_when(
        str_detect(Sample, "_thg_") ~ "Th",
        str_detect(Sample, "_tun_") ~ "Tm",
        str_detect(Sample, "_con_") ~ "DMSO",
        TRUE ~ NA_character_
      ),
      
      idx = as.integer(str_match(Sample, "_(?:thg|tun|con)_(\\d+)\\.raw")[, 2]),
      
      Time = case_when(
        (str_detect(Sample, "_thg_") | str_detect(Sample, "_tun_")) & idx %in% 1:3   ~ "6h",
        (str_detect(Sample, "_thg_") | str_detect(Sample, "_tun_")) & idx %in% 4:6   ~ "12h",
        (str_detect(Sample, "_thg_") | str_detect(Sample, "_tun_")) & idx %in% 7:9   ~ "24h",
        (str_detect(Sample, "_thg_") | str_detect(Sample, "_tun_")) & idx %in% 10:12 ~ "48h",
        
        str_detect(Sample, "_con_") & idx %in% 4:6   ~ "6h",
        str_detect(Sample, "_con_") & idx %in% 7:9   ~ "12h",
        str_detect(Sample, "_con_") & idx %in% 10:12 ~ "24h",
        str_detect(Sample, "_con_") & idx %in% 13:15 ~ "48h",
        
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Metric), !is.na(Treatment), !is.na(Time))
  
  summary_df <- long_df %>%
    group_by(PG.ProteinAccessions, PG.Genes, Treatment, Time, Metric) %>%
    summarise(
      n_rep = sum(!is.na(Intensity)),
      Mean = ifelse(
        n_rep >= min_rep,
        mean(Intensity, na.rm = TRUE),
        NA_real_
      ),
      CV = ifelse(
        n_rep >= min_rep,
        sd(Intensity, na.rm = TRUE) / mean(Intensity, na.rm = TRUE) * 100,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::select(-n_rep) %>%   
    pivot_wider(
      names_from = c(Time, Treatment, Metric),
      values_from = c(Mean, CV),
      names_sep = "_"
    )
  
  summary_df
}

df1_summary <- calculate_mean_cv_4metrics(df1_with_metrics, min_rep = 3)

write_csv(df1_summary, "D:/personal/UVA/SILAC/rebuttal/df1_summary.csv")

df1_summary<-read.csv("D:/personal/UVA/SILAC/rebuttal/df1_summary.csv")

df1_with_metrics_plus_summary <- df1_with_metrics %>%
  left_join(df1_summary, by = c("PG.ProteinAccessions", "PG.Genes"))

write_csv(df1_with_metrics_plus_summary, "D:/personal/UVA/SILAC/rebuttal/df1_with_metrics_plus_summary.csv")

df1_with_metrics_plus_summary<-read.csv("D:/personal/UVA/SILAC/rebuttal/df1_with_metrics_plus_summary.csv")

#Plot one protein as example
# ----
plot_protein_4metrics_bar_stack <- function(df, protein_id = "P11021", min_rep = 3) {
  
  id_cols <- c("PG.ProteinAccessions", "PG.Genes")
  
  raw_cols <- names(df)[grepl(
    "PG\\.MS2Channel1$|PG\\.MS2Channel2$|PG\\.Total$|PG\\.HL_ratio$",
    names(df)
  )]
  
  dat0 <- df %>%
    dplyr::select(dplyr::any_of(id_cols), dplyr::all_of(raw_cols)) %>%
    dplyr::filter(stringr::str_detect(PG.ProteinAccessions, stringr::fixed(protein_id)))
  
  if (nrow(dat0) == 0) stop("Protein not found in df.")
  
  gene_title <- suppressWarnings(stats::na.omit(dat0$PG.Genes)[1])
  if (length(gene_title) == 0 || is.na(gene_title) || gene_title == "") gene_title <- protein_id
  
  long_df <- dat0 %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(id_cols),
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    dplyr::mutate(
      Intensity = suppressWarnings(as.numeric(Intensity)),
      
      Metric = dplyr::case_when(
        stringr::str_detect(Sample, "PG\\.MS2Channel1$") ~ "Light",
        stringr::str_detect(Sample, "PG\\.MS2Channel2$") ~ "Heavy",
        stringr::str_detect(Sample, "PG\\.Total$") ~ "Total",
        stringr::str_detect(Sample, "PG\\.HL_ratio$") ~ "HL",
        TRUE ~ NA_character_
      ),
      
      Treatment = dplyr::case_when(
        stringr::str_detect(Sample, "_thg_") ~ "Th",
        stringr::str_detect(Sample, "_tun_") ~ "Tm",
        stringr::str_detect(Sample, "_con_") ~ "DMSO",
        TRUE ~ NA_character_
      ),
      
      idx = as.integer(stringr::str_match(Sample, "_(?:thg|tun|con)_(\\d+)\\.raw")[, 2]),
      
      Time = dplyr::case_when(
        (stringr::str_detect(Sample, "_thg_") | stringr::str_detect(Sample, "_tun_")) & idx %in% 1:3   ~ "6h",
        (stringr::str_detect(Sample, "_thg_") | stringr::str_detect(Sample, "_tun_")) & idx %in% 4:6   ~ "12h",
        (stringr::str_detect(Sample, "_thg_") | stringr::str_detect(Sample, "_tun_")) & idx %in% 7:9   ~ "24h",
        (stringr::str_detect(Sample, "_thg_") | stringr::str_detect(Sample, "_tun_")) & idx %in% 10:12 ~ "48h",
        
        stringr::str_detect(Sample, "_con_") & idx %in% 4:6   ~ "6h",
        stringr::str_detect(Sample, "_con_") & idx %in% 7:9   ~ "12h",
        stringr::str_detect(Sample, "_con_") & idx %in% 10:12 ~ "24h",
        stringr::str_detect(Sample, "_con_") & idx %in% 13:15 ~ "48h",
        
        TRUE ~ NA_character_
      ),
      
      Replicate = dplyr::case_when(
        idx %in% c(1, 4, 7, 10, 13) ~ 1L,
        idx %in% c(2, 5, 8, 11, 14) ~ 2L,
        idx %in% c(3, 6, 9, 12, 15) ~ 3L,
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::filter(!is.na(Metric), !is.na(Treatment), !is.na(Time), !is.na(Replicate))
  
  # Strict mean/sd/sem: require >= min_rep non-NA replicates
  sum_df <- long_df %>%
    dplyr::group_by(Time, Treatment, Metric) %>%
    dplyr::summarise(
      n = sum(!is.na(Intensity)),
      mean_int = if (n >= min_rep) mean(Intensity, na.rm = TRUE) else NA_real_,
      sd = if (n >= min_rep) sd(Intensity, na.rm = TRUE) else NA_real_,
      sem = if (n >= min_rep) sd(Intensity, na.rm = TRUE) / sqrt(n) else NA_real_,
      .groups = "drop"
    )
  
  p_to_star <- function(p) {
    dplyr::case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ ""
    )
  }
  
  # Strict log2FC and p-value: require >= min_rep non-NA in both treated and control
  test_df <- long_df %>%
    dplyr::group_by(Time, Metric) %>%
    dplyr::group_modify(~{
      dat <- .x
      control_vals <- dat %>% dplyr::filter(Treatment == "DMSO") %>% dplyr::pull(Intensity)
      n_ct <- sum(!is.na(control_vals))
      
      dplyr::bind_rows(lapply(c("Th", "Tm"), function(tr) {
        tr_vals <- dat %>% dplyr::filter(Treatment == tr) %>% dplyr::pull(Intensity)
        n_tr <- sum(!is.na(tr_vals))
        
        mean_tr <- if (n_tr >= min_rep) mean(tr_vals, na.rm = TRUE) else NA_real_
        mean_ct <- if (n_ct >= min_rep) mean(control_vals, na.rm = TRUE) else NA_real_
        
        log2fc <- if (!is.na(mean_tr) && !is.na(mean_ct) && mean_tr > 0 && mean_ct > 0) {
          log2(mean_tr / mean_ct)
        } else {
          NA_real_
        }
        
        pval <- NA_real_
        if (n_tr >= min_rep && n_ct >= min_rep) {
          pval <- tryCatch(stats::t.test(tr_vals, control_vals)$p.value, error = function(e) NA_real_)
        }
        
        tibble::tibble(
          Treatment = tr,
          n_tr = n_tr,
          n_ct = n_ct,
          log2FC = log2fc,
          p_value = pval,
          p_star = p_to_star(pval)
        )
      }))
    }) %>%
    dplyr::ungroup()
  
  plot_df <- sum_df %>%
    dplyr::left_join(test_df, by = c("Time", "Metric", "Treatment")) %>%
    dplyr::mutate(
      Time = factor(Time, levels = c("6h", "12h", "24h", "48h")),
      Treatment = factor(Treatment, levels = c("DMSO", "Tm", "Th")),
      Metric = factor(Metric, levels = c("Light", "Heavy", "Total", "HL")),
      sem = ifelse(is.na(sem), 0, sem),
      fc_num = dplyr::case_when(
        Treatment == "DMSO" ~ "",
        is.na(log2FC) ~ "NA",
        TRUE ~ sprintf("%.2f", log2FC)
      ),
      label_txt = dplyr::case_when(
        Treatment == "DMSO" ~ "",
        TRUE ~ paste0(fc_num, ifelse(p_star == "", "", paste0("\n", p_star)))
      )
    )
  
  plot_df <- plot_df %>%
    dplyr::group_by(Metric) %>%
    dplyr::mutate(
      metric_max = max(mean_int + sem, na.rm = TRUE),
      offset = 0.035 * ifelse(is.finite(metric_max) & metric_max > 0, metric_max, 1),
      y_lab = mean_int + sem + offset
    ) %>%
    dplyr::ungroup()
  
  base_theme <- ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10, face = NULL),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    )
  
  fill_scale <- ggplot2::scale_fill_manual(values = c(
    "DMSO" = "#1F81BD",
    "Tm"   = "#EB8124",
    "Th"   = "#B9B9B9"
  ))
  
  make_one_plot <- function(metric_name, show_title = FALSE) {
    d <- plot_df %>% dplyr::filter(Metric == metric_name)
    
    ymax <- max(d$y_lab, d$mean_int + d$sem, na.rm = TRUE)
    ymax <- ifelse(is.finite(ymax), ymax * 1.08, 1)
    
    ggplot2::ggplot(d, ggplot2::aes(x = Time, y = mean_int, fill = Treatment)) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.85),
        width = 0.75,
        color = NA
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mean_int - sem, ymax = mean_int + sem),
        position = ggplot2::position_dodge(width = 0.8),
        width = 0.15
      ) +
      ggplot2::geom_text(
        data = d %>% dplyr::filter(Treatment != "DMSO"),
        ggplot2::aes(x = Time, y = y_lab, label = label_txt, group = Treatment),
        position = ggplot2::position_dodge(width = 0.8),
        size = 3.5,
        vjust = 0,
        lineheight = 0.9
      ) +
      ggplot2::labs(
        title = if (show_title) gene_title else NULL,
        y = metric_name
      ) +
      fill_scale +
      ggplot2::coord_cartesian(ylim = c(0, ymax), clip = "off") +
      base_theme +
      ggplot2::theme(
        plot.title = if (show_title) ggplot2::element_text(size = 10, hjust = 0.5) else ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 10, 5, 5)
      )
  }
  
  p_light <- make_one_plot("Light", show_title = TRUE)
  p_heavy <- make_one_plot("Heavy")
  p_total <- make_one_plot("Total")
  p_hl    <- make_one_plot("HL")
  
  p_stack <- (p_light | p_heavy) / (p_total | p_hl) +
    patchwork::plot_layout(heights = c(1, 1))
  
  list(
    plot = p_stack,
    plot_df = plot_df,
    stats_vs_dmso = test_df,
    gene_title = gene_title
  )
}

# Run for P11021
res_P11021  <- plot_protein_4metrics_bar_stack(
  df1_with_metrics_plus_summary,
  protein_id = "P11021",
  min_rep = 3
)

# Show stacked plot
res_P11021$plot
res_P11021$stats_vs_dmso
print(res_P11021$stats_vs_dmso, n = Inf)
# ----



#CV boxplot
#----
cond_fill <- c(
  "DMSO" = "#1F81BD", "Tm" = "#EB8124", "Th" = "#B9B9B9"
)
times <- c("6h", "12h", "24h", "48h")
conditions <- c("DMSO", "Tm","Th" )
metrics <- c("Light", "Heavy", "Total", "HL")

# Convert Mean columns to long format
mean_long <- df1_summary %>%
  pivot_longer(
    cols = starts_with("Mean_"),
    names_to = c("Time", "Condition", "Metric"),
    names_pattern = "Mean_(.*)_(.*)_(.*)",
    values_to = "Mean"
  )

# Keep only 4 time points
mean_long <- mean_long %>%
  filter(Time %in% times)

# Keep proteins that have non-NA mean at all 4 time points
proteins_keep <- mean_long %>%
  group_by(PG.ProteinAccessions, Condition, Metric) %>%
  summarise(n_nonNA = sum(!is.na(Mean)), .groups = "drop") %>%
  filter(n_nonNA == 4)


cv_long <- df1_summary %>%
  pivot_longer(
    cols = starts_with("CV_"),
    names_to = c("Time", "Condition", "Metric"),
    names_pattern = "CV_(.*)_(.*)_(.*)",
    values_to = "CV"
  ) %>%
  filter(Time %in% times)

cv_filtered <- cv_long %>%
  inner_join(proteins_keep,
             by = c("PG.ProteinAccessions", "Condition", "Metric"))

plot_cv_metric <- function(metric_name) {
  data_subset <- cv_filtered %>%
    filter(Metric == metric_name, Time %in% times, Condition %in% conditions) %>%
    mutate(
      Time = factor(Time, levels = times),
      Condition = factor(Condition, levels = conditions)
    )
  
  ggplot(data_subset, aes(x = Time, y = CV)) +
    geom_boxplot(
      width = 0.75,
      color = "black",
      aes(fill = Condition),
      outlier.shape = NA
    ) +
    #stat_summary(fun = mean, geom = "point", color = "black", shape = 18, size = 2.6) +
    geom_hline(yintercept = 20, linetype = "dashed", color= "grey",linewidth = 0.6) +
    coord_cartesian(ylim = c(0, 50)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = cond_fill, drop = FALSE) +
    facet_wrap(~ Condition, nrow = 1) +
    labs(title = metric_name, x = "Time", y = "CV (%)") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.line = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.18, "cm"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 10)
    )
}

p_light <- plot_cv_metric("Light")
p_heavy <- plot_cv_metric("Heavy")
p_total <- plot_cv_metric("Total")
p_hl    <- plot_cv_metric("HL")

# Arrange 4 metric plots into 2 rows
p_all <- (p_light | p_heavy) / (p_total | p_hl)
p_all
#----

#protein count bar plot

# Count protein coverage per condition + 3-condition overlap
count_protein_timepoint_coverage_with_overlap <- function(df_summary) {
  
  mean_cols <- names(df_summary)[grepl("^Mean_", names(df_summary))]
  
  long_mean <- df_summary %>%
    dplyr::select(PG.ProteinAccessions, PG.Genes, all_of(mean_cols)) %>%
    pivot_longer(
      cols = all_of(mean_cols),
      names_to = "name",
      values_to = "value"
    ) %>%
    mutate(
      Time      = str_match(name, "^Mean_([^_]+)_([^_]+)_([^_]+)$")[, 2],
      Condition = str_match(name, "^Mean_([^_]+)_([^_]+)_([^_]+)$")[, 3],
      Metric    = str_match(name, "^Mean_([^_]+)_([^_]+)_([^_]+)$")[, 4]
    ) %>%
    filter(!is.na(Time), !is.na(Condition), !is.na(Metric))
  
  # Number of non-NA timepoints per protein within each Condition x Metric
  protein_coverage <- long_mean %>%
    group_by(PG.ProteinAccessions, PG.Genes, Condition, Metric) %>%
    summarise(
      n_timepoints_nonNA = n_distinct(Time[!is.na(value)]),
      .groups = "drop"
    )
  
  # Convert to long thresholds (>=2, >=3, >=4)
  protein_pass_long <- protein_coverage %>%
    tidyr::crossing(threshold = c(2L, 3L, 4L)) %>%
    mutate(pass = n_timepoints_nonNA >= threshold)
  
  # 1) Per-condition counts
  per_condition_counts <- protein_pass_long %>%
    filter(pass) %>%
    group_by(Condition, Metric, threshold) %>%
    summarise(
      n_proteins = n_distinct(PG.ProteinAccessions),
      .groups = "drop"
    ) %>%
    mutate(
      Condition = factor(Condition, levels = c("DMSO", "Th", "Tm")),
      Metric = factor(Metric, levels = c("Light", "Heavy", "Total", "HL"))
    )
  
  # 2) 3-condition overlap counts:
  # protein must pass in all three conditions for the same Metric + threshold
  overlap_counts <- protein_pass_long %>%
    filter(pass, Condition %in% c("DMSO", "Th", "Tm")) %>%
    group_by(PG.ProteinAccessions, Metric, threshold) %>%
    summarise(
      n_conditions_pass = n_distinct(Condition),
      .groups = "drop"
    ) %>%
    filter(n_conditions_pass == 3) %>%
    group_by(Metric, threshold) %>%
    summarise(
      n_proteins = n_distinct(PG.ProteinAccessions),
      .groups = "drop"
    ) %>%
    mutate(
      Condition = "Overlap_3cond",
      Condition = factor(Condition, levels = c("DMSO", "Th", "Tm", "Overlap_3cond")),
      Metric = factor(Metric, levels = c("Light", "Heavy", "Total", "HL"))
    )
  
  # Combine
  combined_long <- bind_rows(
    per_condition_counts,
    overlap_counts
  ) %>%
    arrange(Metric, threshold, Condition)
  
  # Wide table (one row per Metric x threshold, columns = DMSO/Th/Tm/Overlap)
  combined_table <- combined_long %>%
    mutate(threshold_label = paste0(">=", threshold, " timepoints")) %>%
    dplyr::select(Metric, threshold, threshold_label, Condition, n_proteins) %>%
    pivot_wider(
      names_from = Condition,
      values_from = n_proteins
    ) %>%
    arrange(factor(Metric, levels = c("Light", "Heavy", "Total", "HL")), desc(threshold))
  
  # Optional prettier order of threshold rows: 4,3,2
  combined_table <- combined_table %>%
    mutate(threshold = factor(threshold, levels = c(4, 3, 2))) %>%
    arrange(Metric, threshold) %>%
    mutate(threshold = as.integer(as.character(threshold)))
  
  list(
    protein_coverage = protein_coverage,
    per_condition_counts = per_condition_counts,
    overlap_counts = overlap_counts,
    combined_long = combined_long,
    combined_table = combined_table
  )
}

# Run
res_cov2 <- count_protein_timepoint_coverage_with_overlap(df1_summary)

#Before CV20 filtering
# -----
protein_cov <- res_cov2$protein_coverage %>%
  mutate(
    Condition = factor(Condition, levels = c("DMSO", "Th", "Tm")),
    Metric = factor(Metric, levels = c("Light", "Heavy", "Total", "HL"))
  )

threshold_n <- 4L

# 1) Bar plot data (threshold = 4)
# Keep Light / Heavy / Total only
# Total and HL are treated as one combined category in plotting
bar_df_4 <- protein_cov %>%
  filter(Metric %in% c("Light", "Heavy", "Total")) %>%
  group_by(Metric, Condition) %>%
  summarise(
    count = sum(n_timepoints_nonNA >= threshold_n, na.rm = TRUE),
    .groups = "drop"
  )

pretty_metric <- c(
  "Light" = "Light",
  "Heavy" = "Heavy",
  "Total" = "Total (H/L ratio)"
)


make_bar_plot_4 <- function(metric_name) {
  d <- bar_df_4 %>% filter(Metric == metric_name)
  
  ggplot(d, aes(x = Condition, y = count, fill = Condition)) +
    geom_col(width = 0.6, color = NA) +
    geom_text(aes(label = count), vjust = -0.25, size = 3) +
    scale_fill_manual(values = c("DMSO" = "#1F81BD", "Tm" = "#EB8124", "Th" = "#B9B9B9")) +
    labs(
      title = paste0(pretty_metric[[as.character(metric_name)]], " (>=4 timepoints)"),
      x = NULL,
      y = "Protein count"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "none",
      plot.margin = margin(5, 8, 5, 8)
    ) +
    expand_limits(y = max(d$count, na.rm = TRUE) * 1.12)
}

#facet
bar_facet_df_4 <- bar_df_4 %>%
  filter(Metric %in% c("Light", "Heavy", "Total")) %>%
  mutate(
    Condition = factor(Condition, levels = c("DMSO", "Tm", "Th")),
    Metric_label = factor(
      recode(Metric, "Total" = "Total (H/L ratio)"),
      levels = c("Light", "Heavy", "Total (H/L ratio)")
    )
  )

p_bar_facet <- ggplot(bar_facet_df_4, aes(x = Condition, y = count, fill = Condition)) +
  geom_col(width = 0.75, color = NA) +
  geom_text(aes(label = count), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("DMSO" = "#1F81BD", "Tm" = "#EB8124", "Th" = "#B9B9B9")) +
  facet_wrap(~ Metric_label, nrow = 1) +
  labs(
    title = "Protein coverage (>=4 timepoints)",
    x = NULL,
    y = "Protein count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.margin = margin(5, 8, 5, 8)
  ) +
  expand_limits(y = max(bar_facet_df_4$count, na.rm = TRUE) * 1.12)





# 2) Overlap bar plot (3-condition overlap) across metrics
# Keep Light / Heavy / Total only
overlap_df_4 <- protein_cov %>%
  filter(Metric %in% c("Light", "Heavy", "Total")) %>%
  mutate(pass = n_timepoints_nonNA >= threshold_n) %>%
  filter(Condition %in% c("DMSO", "Th", "Tm")) %>%
  group_by(PG.ProteinAccessions, Metric) %>%
  summarise(
    n_cond_pass = sum(pass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Metric) %>%
  summarise(
    overlap_count = sum(n_cond_pass == 3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Light", "Heavy", "Total")),
    Metric_label = factor(
      recode(as.character(Metric),
             "Light" = "Light",
             "Heavy" = "Heavy",
             "Total" = "Total (H/L ratio)"),
      levels = c("Light", "Heavy", "Total (H/L ratio)")
    )
  ) %>%
  arrange(Metric)

p_bar_overlap <- ggplot(overlap_df_4, aes(x = Metric_label, y = overlap_count)) +
  geom_col(width = 0.75, fill = "#6C63B5", color = NA) +
  geom_text(aes(label = overlap_count), vjust = -0.25, size = 3) +
  labs(
    title = "3-condition overlap",
    subtitle = ">=4 timepoints in DMSO, Th, and Tm",
    x = NULL,
    y = "Protein count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 8, 5, 8)
  ) +
  expand_limits(y = max(overlap_df_4$overlap_count, na.rm = TRUE) * 1.12)


p_cov4 <- p_bar_facet | p_bar_overlap
p_cov4 + plot_layout(widths = c(3, 1))


#----

#After CV20% filtering
#----
# 1) Parse Mean and CV columns from df1_summary
mean_long <- df1_summary %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes, starts_with("Mean_")) %>%
  pivot_longer(
    cols = starts_with("Mean_"),
    names_to = "name",
    values_to = "Mean"
  ) %>%
  mutate(
    Time = str_match(name, "^Mean_([^_]+)_([^_]+)_([^_]+)$")[, 2],
    Condition = str_match(name, "^Mean_([^_]+)_([^_]+)_([^_]+)$")[, 3],
    Metric = str_match(name, "^Mean_([^_]+)_([^_]+)_([^_]+)$")[, 4]
  ) %>%
  dplyr::select(-name)

cv_long <- df1_summary %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes, starts_with("CV_")) %>%
  pivot_longer(
    cols = starts_with("CV_"),
    names_to = "name",
    values_to = "CV"
  ) %>%
  mutate(
    Time = str_match(name, "^CV_([^_]+)_([^_]+)_([^_]+)$")[, 2],
    Condition = str_match(name, "^CV_([^_]+)_([^_]+)_([^_]+)$")[, 3],
    Metric = str_match(name, "^CV_([^_]+)_([^_]+)_([^_]+)$")[, 4]
  ) %>%
  dplyr::select(-name)

summary_long <- mean_long %>%
  left_join(
    cv_long,
    by = c("PG.ProteinAccessions", "PG.Genes", "Time", "Condition", "Metric")
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("DMSO", "Th", "Tm")),
    Metric = factor(Metric, levels = c("Light", "Heavy", "Total", "HL")),
    pass = !is.na(Mean) & !is.na(CV) & CV < 20
  )

protein_pass <- summary_long %>%
  group_by(PG.ProteinAccessions, PG.Genes, Condition, Metric) %>%
  summarise(
    n_timepoints_pass = sum(pass, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pass_4tp = n_timepoints_pass == 4)

plot_metrics <- c("Light", "Heavy", "Total", "HL")

bar_df <- protein_pass %>%
  filter(Metric %in% plot_metrics) %>%
  group_by(Metric, Condition) %>%
  summarise(
    count = sum(pass_4tp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = plot_metrics)
  )

pretty_metric <- c(
  "Light" = "Light",
  "Heavy" = "Heavy",
  "Total" = "Total",
  "HL"    = "H/L ratio"
)

bar_df <- bar_df %>%
  mutate(Metric_label = recode(as.character(Metric), !!!pretty_metric))

# 3)Keep 4 metrics separately
plot_metrics <- c("Light", "Heavy", "HL", "Total")

pretty_metric <- c(
  "Light" = "Light",
  "Heavy" = "Heavy",
  "HL"    = "H/L ratio",
  "Total" = "Total"
)

# 1) Facet bar plot for each condition
bar_df_facet <- bar_df %>%
  filter(Metric %in% plot_metrics) %>%
  mutate(
    Condition = factor(Condition, levels = c("DMSO", "Tm", "Th")),
    Metric = factor(Metric, levels = plot_metrics),
    Metric_label = recode(as.character(Metric), !!!pretty_metric),
    Metric_label = factor(
      Metric_label,
      levels = pretty_metric[plot_metrics]
    )
  )

p_bar_facet <- ggplot(bar_df_facet, aes(x = Condition, y = count, fill = Condition)) +
  geom_col(width = 0.75, color = NA) +
  geom_text(aes(label = count), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("DMSO" = "#1F81BD", "Tm" = "#EB8124", "Th" = "#B9B9B9")) +
  facet_wrap(~Metric_label, nrow = 1) +
  labs(
    title = "Protein coverage",
    subtitle = "4 timepoints with non-NA Mean and CV < 20",
    x = NULL,
    y = "Protein count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.margin = margin(5, 8, 5, 8)
  ) +
  expand_limits(y = max(bar_df_facet$count, na.rm = TRUE) * 1.12)

# 2) 3-condition overlap counts
# Protein must pass in DMSO, Th, and Tm for the same metric
overlap_df <- protein_pass %>%
  filter(Metric %in% plot_metrics) %>%
  group_by(PG.ProteinAccessions, Metric) %>%
  summarise(
    n_cond_pass = sum(pass_4tp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Metric) %>%
  summarise(
    overlap_count = sum(n_cond_pass == 3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = plot_metrics),
    Metric_label = recode(as.character(Metric), !!!pretty_metric),
    Metric_label = factor(
      Metric_label,
      levels = pretty_metric[plot_metrics]
    )
  ) %>%
  arrange(Metric)

p_bar_overlap <- ggplot(overlap_df, aes(x = Metric_label, y = overlap_count)) +
  geom_col(width = 0.75, fill = "#6C63B5", color = NA) +
  geom_text(aes(label = overlap_count), vjust = -0.25, size = 3) +
  scale_y_continuous(
    limits = c(0, max(bar_df_facet$count, na.rm = TRUE) * 1.12),
    breaks = c(0, 2000, 4000)
  ) +
  labs(
    title = "overlap",
    subtitle = NULL,
    x = NULL,
    y = "Protein count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8, color = "black", angle = 0, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(2, "mm"),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 8)
  )

# 3) Combine
p_final <- p_bar_facet | p_bar_overlap
p_final + plot_layout(widths = c(4, 1.2))

p<- p_cov4 / p_final
print(p)
#----






