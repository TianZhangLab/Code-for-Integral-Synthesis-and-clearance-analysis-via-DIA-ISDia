library(org.Hs.eg.db)  
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)
library(minpack.lm)
library(purrr)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(plotly)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(patchwork)

#Hela WT combined
df3 <- read.delim("D:/personal/UVA/Data analysis/202412 SILAC/DIA combination/20250402_121853_0327_A012_IT3.5_combine_Report.tsv", header = TRUE, sep = "\t")
df3 <- df3  %>%
  rename_with(
    ~gsub("X\\.[0-9]+\\.\\.[0-9]+_(con|tun|thg)([0-9]).*MS2Channel([0-9])", "\\1\\2_Channel\\3", .),
    starts_with("X.")
  )

#Hela PERK KO combined
df4 <- read.delim("D:/personal/UVA/Data analysis/202412 SILAC/DIA combination/20250420_162456_H730_SILAC_IT3.5_combined_Report.tsv", header = TRUE, sep = "\t")
df4 <- df4  %>%
  rename_with(
    ~gsub("X\\.[0-9]+\\.\\.[0-9]+_(con|tun|thg)([0-9]).*MS2Channel([0-9])", "\\1\\2_Channel\\3", .),
    starts_with("X.")
  )


#calculation
safe_mean3 <- function(x) {
  if (sum(!is.na(x)) == 3) mean(x, na.rm = TRUE) else NA_real_
}
safe_se3 <- function(x) {
  if (sum(!is.na(x)) == 3) sd(x, na.rm = TRUE) / sqrt(3) else NA_real_
}
safe_sum2 <- function(x, y) {
  if_else(!is.na(x) & !is.na(y), x + y, NA_real_)
}
safe_ratio <- function(num, den) {
  if_else(!is.na(num) & !is.na(den), num / den, NA_real_)
}
safe_cv <- function(vec) {
  if (sum(!is.na(vec)) >= 3) {
    sd(vec, na.rm = TRUE) / mean(vec, na.rm = TRUE) * 100} else {
      NA_real_}}
safe_ttest_3 <- function(x, y) {
  if (sum(!is.na(x)) >= 3 && sum(!is.na(y)) >= 3) {
    tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}
calculate_cv <- function(df) {
  df %>%
    mutate(
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
    rowwise() %>%
    mutate(
      # Mean & SE
      con_light_mean = safe_mean3(c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      tun_light_mean = safe_mean3(c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1))),
      thg_light_mean = safe_mean3(c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1))),
      
      con_light_SE = safe_se3(c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      tun_light_SE = safe_se3(c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1))),
      thg_light_SE = safe_se3(c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1))),
      
      con_heavy_mean = safe_mean3(c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      tun_heavy_mean = safe_mean3(c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2))),
      thg_heavy_mean = safe_mean3(c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2))),
      
      con_heavy_SE = safe_se3(c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      tun_heavy_SE = safe_se3(c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2))),
      thg_heavy_SE = safe_se3(c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2))),
      
      con_total_mean = safe_mean3(c_across(c(con1_total, con2_total, con3_total))),
      tun_total_mean = safe_mean3(c_across(c(tun1_total, tun2_total, tun3_total))),
      thg_total_mean = safe_mean3(c_across(c(thg1_total, thg2_total, thg3_total))),
      
      con_total_SE = safe_se3(c_across(c(con1_total, con2_total, con3_total))),
      tun_total_SE = safe_se3(c_across(c(tun1_total, tun2_total, tun3_total))),
      thg_total_SE = safe_se3(c_across(c(thg1_total, thg2_total, thg3_total))),
      
      con_HL_mean = safe_mean3(c_across(c(con1_HL, con2_HL, con3_HL))),
      tun_HL_mean = safe_mean3(c_across(c(tun1_HL, tun2_HL, tun3_HL))),
      thg_HL_mean = safe_mean3(c_across(c(thg1_HL, thg2_HL, thg3_HL))),
      
      con_HL_SE = safe_se3(c_across(c(con1_HL, con2_HL, con3_HL))),
      tun_HL_SE = safe_se3(c_across(c(tun1_HL, tun2_HL, tun3_HL))),
      thg_HL_SE = safe_se3(c_across(c(thg1_HL, thg2_HL, thg3_HL))),
      
      # CVs
      con_light_CV = safe_cv(c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))),
      tun_light_CV = safe_cv(c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1))),
      thg_light_CV = safe_cv(c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1))),
      
      con_heavy_CV = safe_cv(c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))),
      tun_heavy_CV = safe_cv(c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2))),
      thg_heavy_CV = safe_cv(c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2))),
      
      con_total_CV = safe_cv(c_across(c(con1_total, con2_total, con3_total))),
      tun_total_CV = safe_cv(c_across(c(tun1_total, tun2_total, tun3_total))),
      thg_total_CV = safe_cv(c_across(c(thg1_total, thg2_total, thg3_total))),
      
      con_HL_CV = safe_cv(c_across(c(con1_HL, con2_HL, con3_HL))),
      tun_HL_CV = safe_cv(c_across(c(tun1_HL, tun2_HL, tun3_HL))),
      thg_HL_CV = safe_cv(c_across(c(thg1_HL, thg2_HL, thg3_HL))),
      
      # log2 fold change
      light_tun_vs_con = log2(tun_light_mean / con_light_mean),
      light_thg_vs_con = log2(thg_light_mean / con_light_mean),
      heavy_tun_vs_con = log2(tun_heavy_mean / con_heavy_mean),
      heavy_thg_vs_con = log2(thg_heavy_mean / con_heavy_mean),
      total_tun_vs_con = log2(tun_total_mean / con_total_mean),
      total_thg_vs_con = log2(thg_total_mean / con_total_mean),
      HL_tun_vs_con = log2(tun_HL_mean / con_HL_mean),
      HL_thg_vs_con = log2(thg_HL_mean / con_HL_mean),
      
      # p-values
      light_tun_vs_con_p = safe_ttest_3(
        c_across(c(tun1_Channel1, tun2_Channel1, tun3_Channel1)),
        c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))
      ),
      light_thg_vs_con_p = safe_ttest_3(
        c_across(c(thg1_Channel1, thg2_Channel1, thg3_Channel1)),
        c_across(c(con1_Channel1, con2_Channel1, con3_Channel1))
      ),
      heavy_tun_vs_con_p = safe_ttest_3(
        c_across(c(tun1_Channel2, tun2_Channel2, tun3_Channel2)),
        c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))
      ),
      heavy_thg_vs_con_p = safe_ttest_3(
        c_across(c(thg1_Channel2, thg2_Channel2, thg3_Channel2)),
        c_across(c(con1_Channel2, con2_Channel2, con3_Channel2))
      ),
      total_tun_vs_con_p = safe_ttest_3(
        c_across(c(tun1_total, tun2_total, tun3_total)),
        c_across(c(con1_total, con2_total, con3_total))
      ),
      total_thg_vs_con_p = safe_ttest_3(
        c_across(c(thg1_total, thg2_total, thg3_total)),
        c_across(c(con1_total, con2_total, con3_total))
      ),
      HL_tun_vs_con_p = safe_ttest_3(
        c_across(c(tun1_HL, tun2_HL, tun3_HL)),
        c_across(c(con1_HL, con2_HL, con3_HL))
      ),
      HL_thg_vs_con_p = safe_ttest_3(
        c_across(c(thg1_HL, thg2_HL, thg3_HL)),
        c_across(c(con1_HL, con2_HL, con3_HL))
      )
    ) %>%
    ungroup()
}
df3_cv <- calculate_cv(df3)
df4_cv <- calculate_cv(df4)

plot_HL_vs_total_by_group <- function(df, df_name, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  

  light_fc_col  <- paste0("light_",  treatment, "_vs_con")
  light_p_col   <- paste0("light_",  treatment, "_vs_con_p")
  heavy_fc_col  <- paste0("heavy_",  treatment, "_vs_con")
  heavy_p_col   <- paste0("heavy_",  treatment, "_vs_con_p")
  total_fc_col  <- paste0("total_",  treatment, "_vs_con")
  total_p_col   <- paste0("total_",  treatment, "_vs_con_p")
  hl_col        <- paste0("HL_",     treatment, "_vs_con")
  
  # 构建分组并过滤
  df_plot <- df %>%
    select(PG.ProteinAccessions, PG.Genes,
           !!sym(light_fc_col), !!sym(light_p_col),
           !!sym(heavy_fc_col), !!sym(heavy_p_col),
           !!sym(total_fc_col), !!sym(total_p_col),
           !!sym(hl_col)) %>%
    rename(
      Light_FC = !!sym(light_fc_col),
      Light_p  = !!sym(light_p_col),
      Heavy_FC = !!sym(heavy_fc_col),
      Heavy_p  = !!sym(heavy_p_col),
      Total_FC = !!sym(total_fc_col),
      Total_p  = !!sym(total_p_col),
      HL_FC    = !!sym(hl_col)
    ) %>%
    mutate(
      ColorGroup = case_when(
        Light_p > 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "D",
        Light_p > 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "G",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "F",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "B",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "C",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "H",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "A",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "E",
        TRUE ~ "NS"
      )
    ) %>%
    filter(!is.na(Light_FC), !is.na(Heavy_FC), !is.na(HL_FC), !is.na(Total_FC))
  
  # 定义颜色映射
  color_map <- c(
    A = "#FF1F5B",  # LightUp + HeavyUp + TotalUp
    B = "#00CD6C",  # LightUp + TotalUp, HeavyNS
    C = "#F28522",  # LightUp + HeavyDown, TotalNS
    D = "#AF58BA",  # HeavyDown + TotalDown, LightNS
    E = "#FFC61E",  # LightDown + HeavyDown + TotalDown
    F = "#8c564b",  # LightDown + TotalDown, HeavyNS
    G = "#009ADE",  # HeavyUp + TotalUp, LightNS
    H = "#FF69B4",  # LightDown + HeavyUp, TotalNS
    NS = "grey80"   # Others
  )
  df_name_map <- list(
    "df3_cv" = "WT",
    "df4_cv" = expression("PERK"^"-/-")
  )
  
  treat_label <- ifelse(treatment == "tun", "Tm", "Th")
  
  title_expr <- if (df_name %in% names(df_name_map)) {
    substitute(group ~ x * "/DMSO", list(
      group = df_name_map[[df_name]],
      x = treat_label
    ))
  } else {
    paste0(df_name, " ", treat_label, "/DMSO")
  }
  
  cor_res <- cor.test(df_plot$Total_FC, df_plot$HL_FC, use = "complete.obs")
  r_val <- round(cor_res$estimate, 2)
  p_val_display <- ifelse(cor_res$p.value < 1e-16, "< 2.2e-16", formatC(cor_res$p.value, format = "e", digits = 2))
  

  ggplot(df_plot, aes(x = Total_FC, y = HL_FC, color = ColorGroup)) +
    geom_point(
      data  = df_plot %>% filter(ColorGroup == "NS"),
      color = "grey80",
      alpha = 0.5,
      size  = 1.8
    ) +
    geom_point(
      data = df_plot %>% filter(ColorGroup != "NS"),
      aes(color = ColorGroup),
      alpha = 0.8,
      size  = 1.8
    ) +
    geom_hline(yintercept = c(0, 0),
               linetype   = "dashed",
               color      = "grey50",
               linewidth  = 0.8) +
    geom_vline(xintercept = c(-0.58, 0.58),
               linetype   = "dashed",
               color      = "grey50",
               linewidth  = 0.8) +
    annotate("text", x = -3.5, y = 2,
             label = paste0("R = ", r_val, "\nP = ", p_val_display),
             hjust = 0, size = 4, fontface = "bold") +
    scale_color_manual(values = color_map) +
    scale_x_continuous(limits = c(-4, 4)) +
    scale_y_continuous(limits = c(-4, 4)) +
    annotate("text",
             x     = -2.8,
             y     = 3.8,
             label = paste0("n = ", nrow(df_plot)),
             hjust = 1,
             vjust = 1,
             size  = 4) +
    ggtitle(title_expr) +
    labs(
      x     = expression("Total log"[2]*"(FC)"),
      y     = expression("H/L ratio log"[2]*"(FC)"),
      color = NULL
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +  
    theme_minimal()+
    theme(
      panel.grid     = element_blank(),
      panel.border   = element_rect(color = "black", fill = NA),
      axis.line      = element_line(color = "black"),
      axis.text      = element_text(size = 12),
      axis.title     = element_text(size = 14),
      legend.position= "top",
      legend.direction = "horizontal",
      plot.title     = element_text(face = "bold", size = 14, hjust = 0.5)
    )
}
plot_HL_vs_total_by_group(df3_cv,df_name = "df3_cv" ,treatment = "thg")

#heavy vs light FC comparison 
#all colors
plot_light_heavy_WT <- function(df, df_name, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  
  light_fc_col  <- paste0("light_",  treatment, "_vs_con")
  light_p_col   <- paste0("light_",  treatment, "_vs_con_p")
  heavy_fc_col  <- paste0("heavy_",  treatment, "_vs_con")
  heavy_p_col   <- paste0("heavy_",  treatment, "_vs_con_p")
  total_fc_col  <- paste0("total_",  treatment, "_vs_con")
  total_p_col   <- paste0("total_",  treatment, "_vs_con_p")
  hl_col        <- paste0("HL_",     treatment, "_vs_con")

  df_plot <- df %>%
    select(PG.ProteinAccessions, PG.Genes,
           !!sym(light_fc_col), !!sym(light_p_col),
           !!sym(heavy_fc_col), !!sym(heavy_p_col),
           !!sym(total_fc_col), !!sym(total_p_col),
           !!sym(hl_col)) %>%
    rename(
      Light_FC = !!sym(light_fc_col),
      Light_p  = !!sym(light_p_col),
      Heavy_FC = !!sym(heavy_fc_col),
      Heavy_p  = !!sym(heavy_p_col),
      Total_FC = !!sym(total_fc_col),
      Total_p  = !!sym(total_p_col),
      HL_FC    = !!sym(hl_col)
    ) %>%
    mutate(
      ColorGroup = case_when(
        Light_p > 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "D",
        Light_p > 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "G",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "F",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_p > 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "B",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "C",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_p > 0.05 ~ "H",
        Light_FC >  0.58 & Light_p < 0.05 & Heavy_FC >  0.58 & Heavy_p < 0.05 & Total_FC >  0.58 & Total_p < 0.05 ~ "A",
        Light_FC < -0.58 & Light_p < 0.05 & Heavy_FC < -0.58 & Heavy_p < 0.05 & Total_FC < -0.58 & Total_p < 0.05 ~ "E",
        TRUE ~ "NS"
      )
    ) %>%
    filter(!is.na(Light_FC), !is.na(Heavy_FC), !is.na(HL_FC), !is.na(Total_FC))
  
  # 定义颜色映射
  color_map <- c(
    A = "#FF1F5B",  # LightUp + HeavyUp + TotalUp
    B = "#00CD6C",  # LightUp + TotalUp, HeavyNS
    C = "#F28522",  # LightUp + HeavyDown, TotalNS
    D = "#AF58BA",  # HeavyDown + TotalDown, LightNS
    E = "#FFC61E",  # LightDown + HeavyDown + TotalDown
    F = "#8c564b",  # LightDown + TotalDown, HeavyNS
    G = "#009ADE",  # HeavyUp + TotalUp, LightNS
    H = "#FF69B4",  # LightDown + HeavyUp, TotalNS
    NS = "grey80"   # Others
  )
  
  df_name_map <- list(
    "df3_cv" = "WT",
    "df4_cv" = expression("PERK"^"-/-")
  )
  
  treat_label <- ifelse(treatment == "tun", "Tm", "Th")

  title_expr <- if (df_name %in% names(df_name_map)) {
    substitute(group ~ x * "/DMSO", list(
      group = df_name_map[[df_name]],
      x = treat_label
    ))
  } else {
    paste0(df_name, " ", treat_label, "/DMSO")
  }

  cor_res <- cor.test(df_plot$Heavy_FC, df_plot$Light_FC, use = "complete.obs")
  r_val <- round(cor_res$estimate, 2)
  p_val_display <- ifelse(cor_res$p.value < 1e-16, "< 2.2e-16", formatC(cor_res$p.value, format = "e", digits = 2))
  
  ggplot(df_plot, aes(x = Heavy_FC, y = Light_FC, color = ColorGroup)) +
    geom_point(
      data  = df_plot %>% filter(ColorGroup == "NS"),
      color = "grey80",
      alpha = 0.5,
      size  = 1.8
    ) +

    geom_point(
      data = df_plot %>% filter(ColorGroup != "NS"),
      aes(color = ColorGroup),
      alpha = 0.8,
      size  = 1.8
    ) +
    geom_hline(yintercept = c(-0.58, 0.58),
               linetype   = "dashed",
               color      = "grey50",
               linewidth  = 0.8) +
    geom_vline(xintercept = c(-0.58, 0.58),
               linetype   = "dashed",
               color      = "grey50",
               linewidth  = 0.8) +
    annotate("text", x = -3.5, y = 2,
             label = paste0("R = ", r_val, "\nP = ", p_val_display),
             hjust = 0, size = 4, fontface = "bold") +
    scale_color_manual(values = color_map) +
    scale_x_continuous(limits = c(-4, 4)) +
    scale_y_continuous(limits = c(-4, 4)) +
    annotate("text",
             x     = -2.8,
             y     = 3.8,
             label = paste0("n = ", nrow(df_plot)),
             hjust = 1,
             vjust = 1,
             size  = 4) +
    ggtitle(title_expr) +
    labs(
      x = expression("Heavy log"[2]*"(FC)"),
      y = expression("Light log"[2]*"(FC)"),
      color = NULL
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +  
    theme_minimal()+
    theme(
      panel.grid     = element_blank(),
      panel.border   = element_rect(color = "black", fill = NA),
      axis.line      = element_line(color = "black"),
      axis.text      = element_text(size = 12),
      axis.title     = element_text(size = 14),
      legend.position= "top",
      legend.direction = "horizontal",
      plot.title     = element_text(face = "bold", size = 14, hjust = 0.5)
    )
}
plot_light_heavy_WT (df3_cv,df_name = "df3_cv", treatment = "thg")

#pink with H/L>0
df1<-df3_cv%>%
  filter(light_thg_vs_con > 0.58 & light_thg_vs_con_p < 0.05 &
           heavy_thg_vs_con > 0.58 & heavy_thg_vs_con_p  < 0.05 &
           total_thg_vs_con  > 0.58 & total_thg_vs_con_p < 0.05&HL_thg_vs_con>0)

run_go_analysis <- function(df, title_text = "Input") {
  prot_ids <- df$PG.ProteinAccessions
  prot_ids_clean <- sub(";.*", "", prot_ids)
  
  id_map <- bitr(prot_ids_clean, fromType = "UNIPROT",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(id_map) == 0) return(NULL)
  
  go_result <- enrichGO(
    gene = id_map$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  go_filtered <- go_result@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(Source = title_text)
  
  if (nrow(go_filtered) == 0) return(NULL)
  
  return(go_filtered)
}
df1_go <- run_go_analysis(df1)

plot_go_bubble <- function(go_df, top_n = 10, title_text = "GO Enrichment Bubble Plot") {
  if (is.null(go_df)) {
    message("No GO terms to plot.")
    return(NULL)
  }
  
  go_top <- go_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      log_p = -log10(p.adjust),
      Description = factor(Description, levels = rev(Description))  
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  
    scale_fill_gradient(low = "#ffcccc", high = "red", name = "-log10(p.adjust)") +  
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment", 
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(4, 20)) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),  
      axis.line = element_line(color = "black"), 
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
}
plot_go_bubble(df1_go, top_n = 10, title_text = "Top GO Terms (BP)")

#pink with H/L<0
df2<-df3_cv%>%
  filter(light_thg_vs_con > 0.58 & light_thg_vs_con_p < 0.05 &
           heavy_thg_vs_con > 0.58 & heavy_thg_vs_con_p  < 0.05 &
           total_thg_vs_con  > 0.58 & total_thg_vs_con_p < 0.05&HL_thg_vs_con<0)

df2_go <- run_go_analysis(df2)
plot_go_bubble <- function(go_df, top_n = 10, title_text = "GO Enrichment Bubble Plot") {
  if (is.null(go_df)) {
    message("No GO terms to plot.")
    return(NULL)
  }
  
  go_top <- go_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      log_p = -log10(p.adjust),
      Description = factor(Description, levels = rev(Description))  
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  
    scale_fill_gradient(low = "#ffcccc", high = "red", name = "-log10(p.adjust)") +  # 红色渐变
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment",  
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(4, 12)) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),  
      axis.line = element_line(color = "black"), 
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
}
plot_go_bubble(df2_go, top_n = 10, title_text = "Top GO Terms (BP)")


#purple and yellow proteins
df5<-df3_cv%>%
  filter(
    heavy_thg_vs_con < -0.58 & heavy_thg_vs_con_p < 0.05 &
      total_thg_vs_con < -0.58 & total_thg_vs_con_p < 0.05 &
      (
        (light_thg_vs_con_p > 0.05) |
          (light_thg_vs_con < -0.58 & light_thg_vs_con_p < 0.05)
      )
  )

df5_go <- run_go_analysis(df5)
plot_go_bubble <- function(go_df, top_n = 10, title_text = "GO Enrichment Bubble Plot") {
  if (is.null(go_df)) {
    message("No GO terms to plot.")
    return(NULL)
  }
  
  go_top <- go_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      log_p = -log10(p.adjust),
      Description = factor(Description, levels = rev(Description))  
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  
    scale_fill_gradient(low = "blue", high = "blue", name = "-log10(p.adjust)") +  
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment",  
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(0, 11)) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),  
      axis.line = element_line(color = "black"), 
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
}
plot_go_bubble(df5_go, top_n = 10, title_text = "Top GO Terms (BP)")

#uORF protein analysis
uORF_1 <- read_excel("D:/personal/UVA/Data analysis/202412 SILAC/Fig6/uORF_1.xlsx")
uORF_df3<-df3_cv%>%
  filter(PG.ProteinAccessions%in%uORF_1$Entry)
median_val <- median(uORF_df3$heavy_thg_vs_con, na.rm = TRUE)
n_val <- sum(!is.na(uORF_df3$heavy_thg_vs_con))

ggplot(uORF_df3, aes(x = heavy_thg_vs_con)) +
  geom_histogram(binwidth = 0.2, fill = "#b8b8b8", color = "white", boundary = 0) +
  
  geom_vline(xintercept = median_val, linetype = "dashed", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 1.2) +
  scale_y_continuous(limits = c(0,45))+
  scale_x_continuous(limits = c(-2,2))+
  annotate("text",
           x = Inf,
           y = 40,
           label = paste0("Median = ", round(median_val, 2)),
           hjust = 1.1, vjust = 1.5,
           color = "red", size = 4
  ) +
  
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("n = ", n_val),
           hjust = 1.1, vjust = 1.5,
           size = 4
  ) +
  
  labs(
    title = "48h combine",
    x = "Heavy log2(Th/DMSO)",
    y = "Protein count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "none"  
  )

#6h data and uORF
df1 <- read.delim("D:/personal/UVA/Data analysis/202412 SILAC/A012_treated/20250204_132507_A012_SILAC_ALL2_Report with gene.tsv", header = TRUE, sep = "\t")
filter_6h <- function(df) {
  ms2_cols <- grep("PG\\.MS2Channel[12]", colnames(df), value = TRUE)

  ms2_6h <- ms2_cols[
    grepl("thg_[123]\\.raw|tun_[123]\\.raw|con_[456]\\.raw", ms2_cols)
  ]

  df6h <- df %>%
    select(PG.ProteinAccessions, PG.Genes, all_of(ms2_6h))
  
  return(df6h)
}

df1_6h <- filter_6h(df1)

df1_6h <- df1_6h %>%
  rename_with(
    ~ str_replace_all(
      gsub("X\\.[0-9]+\\.\\.[0-9]+_(con|tun|thg)_([0-9]+).*MS2Channel([0-9])", 
           "\\1\\2_Channel\\3", .),
      c("thg1" = "thg1", "thg2" = "thg2", "thg3" = "thg3",
        "tun1" = "tun1", "tun2" = "tun2", "tun3" = "tun3",
        "con4" = "con1", "con5" = "con2", "con6" = "con3")
    )
  )

calculate_FC <- function(df) {
  df %>%
    
    rowwise() %>%
    mutate(
    
      con_heavy_mean = safe_mean3(c(con1_Channel2, con2_Channel2, con3_Channel2)),
      tun_heavy_mean = safe_mean3(c(tun1_Channel2, tun2_Channel2, tun3_Channel2)),
      thg_heavy_mean = safe_mean3(c(thg1_Channel2, thg2_Channel2, thg3_Channel2)),
   
      # log2 fold change
     
      heavy_tun_vs_con  = log2(safe_ratio(tun_heavy_mean, con_heavy_mean)),
      heavy_thg_vs_con  = log2(safe_ratio(thg_heavy_mean, con_heavy_mean)),
      
      # p-values using safe_ttest_3
     
      heavy_tun_vs_con_p  = safe_ttest_3(c(tun1_Channel2, tun2_Channel2, tun3_Channel2),
                                         c(con1_Channel2, con2_Channel2, con3_Channel2)),
      heavy_thg_vs_con_p  = safe_ttest_3(c(thg1_Channel2, thg2_Channel2, thg3_Channel2),
                                         c(con1_Channel2, con2_Channel2, con3_Channel2)),

    ) %>%
    ungroup()
}

df1_6h  <- calculate_FC(df1_6h)

uORF_6h<-df1_6h%>%
  filter(PG.ProteinAccessions%in%uORF_1$Entry)

median_val_6h <- median(uORF_6h$heavy_thg_vs_con, na.rm = TRUE)
n_val_6h <- sum(!is.na(uORF_6h$heavy_thg_vs_con))


ggplot(uORF_6h, aes(x = heavy_thg_vs_con)) +
  geom_histogram(binwidth = 0.2, fill = "#b8b8b8", color = "white", boundary = 0) +
  
  geom_vline(xintercept = median_val, linetype = "dashed", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 1.2) +
  scale_y_continuous(limits = c(0,45))+
  scale_x_continuous(limits = c(-2,2))+
  annotate("text",
           x = Inf,
           y = 40,
           label = paste0("Median = ", round(median_val_6h, 2)),
           hjust = 1.1, vjust = 1.5,
           color = "red", size = 4
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("n = ", n_val_6h),
           hjust = 1.1, vjust = 1.5,
           size = 4
  ) +
  
  labs(
    title = "6h alone",
    x = "Heavy log2(Th/DMSO)",
    y = "Protein count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    legend.position = "none" 
  )

#light_up, heavy_down and total no change
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
    filter(
      .data[[light_col]] > 0.58, .data[[light_p_col]] < 0.05,
      .data[[heavy_col]] < -0.58, .data[[heavy_p_col]] < 0.05,
      .data[[total_p_col]] > 0.05
    ) %>%
    select(PG.ProteinAccessions, PG.Genes)
  
  # df3 subset
  df3_selected <- df3_cv %>%
    filter(PG.ProteinAccessions %in% selected_proteins$PG.ProteinAccessions) %>%
    select(PG.ProteinAccessions, PG.Genes,
           !!light_col, !!heavy_col, !!total_col,
           !!light_p_col, !!heavy_p_col, !!total_p_col)
  
  # df4 subset
  df4_selected <- df4_cv %>%
    filter(PG.ProteinAccessions %in% selected_proteins$PG.ProteinAccessions) %>%
    select(PG.ProteinAccessions, PG.Genes,
           !!light_col, !!heavy_col, !!total_col,
           !!light_p_col, !!heavy_p_col, !!total_p_col)
  
  # rename for merging
  df3_selected <- df3_selected %>%
    rename(
      light_df3 = !!light_col, heavy_df3 = !!heavy_col, total_df3 = !!total_col,
      light_p_df3 = !!light_p_col, heavy_p_df3 = !!heavy_p_col, total_p_df3 = !!total_p_col
    )
  
  df4_selected <- df4_selected %>%
    rename(
      light_df4 = !!light_col, heavy_df4 = !!heavy_col, total_df4 = !!total_col,
      light_p_df4 = !!light_p_col, heavy_p_df4 = !!heavy_p_col, total_p_df4 = !!total_p_col
    )
  
  # merge
  merged_selected <- inner_join(df3_selected, df4_selected, by = c("PG.ProteinAccessions", "PG.Genes"))
  
  # Color group
  merged_selected <- merged_selected %>%
    mutate(ColorGroup = case_when(
      .data[[paste0(channel, "_p_df3")]] < 0.05 & .data[[paste0(channel, "_p_df4")]] >= 0.05 ~ "sig_df3",
      .data[[paste0(channel, "_p_df3")]] >= 0.05 & .data[[paste0(channel, "_p_df4")]] < 0.05 ~ "sig_df4",
      .data[[paste0(channel, "_p_df3")]] < 0.05 & .data[[paste0(channel, "_p_df4")]] < 0.05 ~ "both",
      TRUE ~ "ns"
    ))
  
  # Color map
  color_map <- c("sig_df3" = "#009ade", "sig_df4" = "#ff1f5b", "both" = "#ffc61e", "ns" = "grey80")
  
  # plotting columns
  fc_df3_col <- paste0(channel, "_df3")
  fc_df4_col <- paste0(channel, "_df4")
  
  p <- ggplot(merged_selected, aes(x = .data[[fc_df3_col]], y = .data[[fc_df4_col]], color = ColorGroup)) +
    geom_point(size = 2, alpha = 0.8) +
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
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}
plot_lightup_heavydown_totalns(df3_cv, df4_cv, channel = "light", treatment = "thg")
plot_lightup_heavydown_totalns(df3_cv, df4_cv, channel = "heavy", treatment = "thg")
plot_lightup_heavydown_totalns(df3_cv, df4_cv, channel = "total", treatment = "thg")
