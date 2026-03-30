library(purrr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(plotly)
library(patchwork)

df3 <- read.delim("D:/personal/UVA/Data analysis/202412 SILAC/DIA combination/20250402_121853_0327_A012_IT3.5_combine_Report.tsv", header = TRUE, sep = "\t")
df3 <- df3  %>%
  rename_with(
    ~gsub("X\\.[0-9]+\\.\\.[0-9]+_(con|tun|thg)([0-9]).*MS2Channel([0-9])", "\\1\\2_Channel\\3", .),
    starts_with("X.")
  )

#df3 CV calculate and count
# ---- Start of code block ----
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

df3_cv<-read.csv("D:/personal/UVA/SILAC/rebuttal/df3_cv.csv")

treat_colors<-c("DMSO" = "#1F81BD",
                "Tm"   = "#EB8124",
                "Th"   = "#B9B9B9")

#count protein with 3 replicates
count_proteins_with_3_reps <- function(df) {
  replicate_groups <- list(
    con_light = c("con1_Channel1", "con2_Channel1", "con3_Channel1"),
    con_heavy = c("con1_Channel2", "con2_Channel2", "con3_Channel2"),
    con_total = c("con1_total", "con2_total", "con3_total"),
    con_HL = c("con1_HL", "con2_HL", "con3_HL"),
    
    tun_light = c("tun1_Channel1", "tun2_Channel1", "tun3_Channel1"),
    tun_heavy = c("tun1_Channel2", "tun2_Channel2", "tun3_Channel2"),
    tun_total = c("tun1_total", "tun2_total", "tun3_total"),
    tun_HL = c("tun1_HL", "tun2_HL", "tun3_HL"),
    
    thg_light = c("thg1_Channel1", "thg2_Channel1", "thg3_Channel1"),
    thg_heavy = c("thg1_Channel2", "thg2_Channel2", "thg3_Channel2"),
    thg_total = c("thg1_total", "thg2_total", "thg3_total"),
    thg_HL = c("thg1_HL", "thg2_HL", "thg3_HL")
  )
  
  result <- data.frame(
    Condition = names(replicate_groups),
    Protein_Count = sapply(replicate_groups, function(cols) {
      sum(rowSums(!is.na(df[cols])) == 3)
    })
  )
  
  return(result)
}
cv_summary_df3 <- count_proteins_with_3_reps(df3_cv)

cv_summary_df3 <- cv_summary_df3 %>%
  mutate(
    Treatment = case_when(
      str_detect(Condition, "^con") ~ "DMSO",
      str_detect(Condition, "^tun") ~ "Tm",
      str_detect(Condition, "^thg") ~ "Th",
      TRUE ~ NA_character_
    ),
    Intensity = case_when(
      str_detect(Condition, "light") ~ "Light",
      str_detect(Condition, "heavy") ~ "Heavy",
      str_detect(Condition, "total") ~ "Total(H/L ratio)",
      str_detect(Condition, "HL") ~ "H/L",
      TRUE ~ NA_character_
    )
  )


df_facet <- cv_summary_df3 %>%
  filter(Intensity %in% c("Light", "Heavy", "Total(H/L ratio)")) %>%
  mutate(Intensity = recode(Intensity,
                            "Total(H/L ratio)" = "H/L ratio")) %>%
  mutate(Intensity = factor(Intensity, levels = c("Light", "Heavy", "H/L ratio")))

p_before <-ggplot(df_facet, aes(x = Treatment, y = Protein_Count, fill = Treatment)) +
  geom_col(width = 0.75)+
  geom_text(aes(label = Protein_Count), vjust = -0.3, size = 3, color = "black") +
  scale_fill_manual(values = treat_colors) +
  scale_y_continuous(limits = c(0, 10500),breaks=c(0,5000,10000)) +
  scale_x_discrete(limits = c("DMSO", "Tm", "Th"))+
  facet_wrap(~Intensity, nrow = 1, strip.position = "top") +
  labs(
    x = "Treatment",
    y = "Protein Count"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    theme(panel.spacing = unit(1.5, "cm"))
  )
# overlap counts (3 reps) for each metric
overlap_before <- tibble(
  Intensity = c("Light", "Heavy", "H/L ratio"),
  Protein_Count = c(
    sum(rowSums(!is.na(df3_cv[c("con1_Channel1","con2_Channel1","con3_Channel1")])) == 3 &
          rowSums(!is.na(df3_cv[c("tun1_Channel1","tun2_Channel1","tun3_Channel1")])) == 3 &
          rowSums(!is.na(df3_cv[c("thg1_Channel1","thg2_Channel1","thg3_Channel1")])) == 3),
    sum(rowSums(!is.na(df3_cv[c("con1_Channel2","con2_Channel2","con3_Channel2")])) == 3 &
          rowSums(!is.na(df3_cv[c("tun1_Channel2","tun2_Channel2","tun3_Channel2")])) == 3 &
          rowSums(!is.na(df3_cv[c("thg1_Channel2","thg2_Channel2","thg3_Channel2")])) == 3),
    sum(rowSums(!is.na(df3_cv[c("con1_total","con2_total","con3_total")])) == 3 &
          rowSums(!is.na(df3_cv[c("tun1_total","tun2_total","tun3_total")])) == 3 &
          rowSums(!is.na(df3_cv[c("thg1_total","thg2_total","thg3_total")])) == 3)
  )
) %>%
  mutate(Intensity = factor(Intensity, levels = c("Light","Heavy","H/L ratio")))

p_overlap_before <- ggplot(overlap_before, aes(x = Intensity, y = Protein_Count)) +
  geom_col(width = 0.75, fill = "#6C63B5") +
  geom_text(aes(label = Protein_Count), vjust = -0.3, size = 2.5, color = "black") +
  scale_y_continuous(limits = c(0, 10500), breaks = c(0, 5000, 10000)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(5, 5, 5, 5)
  )

# put overlap panel on the right
p_before + p_overlap_before + plot_layout(widths = c(3, 1))

#count protein in 3 replicates with CV<20%
count_proteins_3reps_cv20 <- function(df, cv_cutoff = 20) {
  groups <- tibble::tribble(
    ~Condition,     ~rep_cols,                                                       ~cv_col,
    "con_light",    c("con1_Channel1","con2_Channel1","con3_Channel1"),              "con_light_CV",
    "con_heavy",    c("con1_Channel2","con2_Channel2","con3_Channel2"),              "con_heavy_CV",
    "con_total",    c("con1_total","con2_total","con3_total"),                       "con_total_CV",
    "con_HL",       c("con1_HL","con2_HL","con3_HL"),                                "con_HL_CV",
    
    "tun_light",    c("tun1_Channel1","tun2_Channel1","tun3_Channel1"),              "tun_light_CV",
    "tun_heavy",    c("tun1_Channel2","tun2_Channel2","tun3_Channel2"),              "tun_heavy_CV",
    "tun_total",    c("tun1_total","tun2_total","tun3_total"),                       "tun_total_CV",
    "tun_HL",       c("tun1_HL","tun2_HL","tun3_HL"),                                "tun_HL_CV",
    
    "thg_light",    c("thg1_Channel1","thg2_Channel1","thg3_Channel1"),              "thg_light_CV",
    "thg_heavy",    c("thg1_Channel2","thg2_Channel2","thg3_Channel2"),              "thg_heavy_CV",
    "thg_total",    c("thg1_total","thg2_total","thg3_total"),                       "thg_total_CV",
    "thg_HL",       c("thg1_HL","thg2_HL","thg3_HL"),                                "thg_HL_CV"
  )
  
  out <- groups %>%
    rowwise() %>%
    mutate(
      Protein_Count = {
        reps_ok <- rowSums(!is.na(df[rep_cols])) == 3
        cv_ok <- !is.na(df[[cv_col]]) & df[[cv_col]] < cv_cutoff
        sum(reps_ok & cv_ok, na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    dplyr::select(Condition, Protein_Count)
  
  out
}

cv_summary_df3 <- count_proteins_3reps_cv20(df3_cv, cv_cutoff = 20) %>%
  mutate(
    Treatment = case_when(
      str_detect(Condition, "^con") ~ "DMSO",
      str_detect(Condition, "^tun") ~ "Tm",
      str_detect(Condition, "^thg") ~ "Th",
      TRUE ~ NA_character_
    ),
    Intensity = case_when(
      str_detect(Condition, "light") ~ "Light",
      str_detect(Condition, "heavy") ~ "Heavy",
      str_detect(Condition, "total") ~ "Total",
      str_detect(Condition, "HL") ~ "H/L ratio",
      TRUE ~ NA_character_
    )
  )

df_facet <- cv_summary_df3 %>%
  filter(Intensity %in% c("Light", "Heavy", "H/L ratio", "Total")) %>%
  mutate(
    Intensity = factor(Intensity, levels = c("Light", "Heavy", "H/L ratio", "Total")),
    Treatment = factor(Treatment, levels = c("DMSO", "Tm", "Th"))
  )

p_after <- ggplot(df_facet, aes(x = Treatment, y = Protein_Count, fill = Treatment)) +
  geom_col(width = 0.75, color = NA) +
  geom_text(aes(label = Protein_Count), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("DMSO" = "#1F81BD", "Tm" = "#EB8124", "Th" = "#B9B9B9")) +
  facet_wrap(~Intensity, nrow = 1) +
  labs(
    title = "Protein coverage",
    subtitle = "3 replicates with CV < 20",
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
  expand_limits(y = max(df_facet$Protein_Count, na.rm = TRUE) * 1.12)


cv_cutoff <- 20

overlap_after <- tibble(
  Intensity = c("Light", "Heavy", "H/L ratio", "Total"),
  Protein_Count = c(
    sum(
      rowSums(!is.na(df3_cv[c("con1_Channel1","con2_Channel1","con3_Channel1")])) == 3 &
        rowSums(!is.na(df3_cv[c("tun1_Channel1","tun2_Channel1","tun3_Channel1")])) == 3 &
        rowSums(!is.na(df3_cv[c("thg1_Channel1","thg2_Channel1","thg3_Channel1")])) == 3 &
        !is.na(df3_cv$con_light_CV) & df3_cv$con_light_CV < cv_cutoff &
        !is.na(df3_cv$tun_light_CV) & df3_cv$tun_light_CV < cv_cutoff &
        !is.na(df3_cv$thg_light_CV) & df3_cv$thg_light_CV < cv_cutoff,
      na.rm = TRUE
    ),
    sum(
      rowSums(!is.na(df3_cv[c("con1_Channel2","con2_Channel2","con3_Channel2")])) == 3 &
        rowSums(!is.na(df3_cv[c("tun1_Channel2","tun2_Channel2","tun3_Channel2")])) == 3 &
        rowSums(!is.na(df3_cv[c("thg1_Channel2","thg2_Channel2","thg3_Channel2")])) == 3 &
        !is.na(df3_cv$con_heavy_CV) & df3_cv$con_heavy_CV < cv_cutoff &
        !is.na(df3_cv$tun_heavy_CV) & df3_cv$tun_heavy_CV < cv_cutoff &
        !is.na(df3_cv$thg_heavy_CV) & df3_cv$thg_heavy_CV < cv_cutoff,
      na.rm = TRUE
    ),
    sum(
      rowSums(!is.na(df3_cv[c("con1_HL","con2_HL","con3_HL")])) == 3 &
        rowSums(!is.na(df3_cv[c("tun1_HL","tun2_HL","tun3_HL")])) == 3 &
        rowSums(!is.na(df3_cv[c("thg1_HL","thg2_HL","thg3_HL")])) == 3 &
        !is.na(df3_cv$con_HL_CV) & df3_cv$con_HL_CV < cv_cutoff &
        !is.na(df3_cv$tun_HL_CV) & df3_cv$tun_HL_CV < cv_cutoff &
        !is.na(df3_cv$thg_HL_CV) & df3_cv$thg_HL_CV < cv_cutoff,
      na.rm = TRUE
    ),
    sum(
      rowSums(!is.na(df3_cv[c("con1_total","con2_total","con3_total")])) == 3 &
        rowSums(!is.na(df3_cv[c("tun1_total","tun2_total","tun3_total")])) == 3 &
        rowSums(!is.na(df3_cv[c("thg1_total","thg2_total","thg3_total")])) == 3 &
        !is.na(df3_cv$con_total_CV) & df3_cv$con_total_CV < cv_cutoff &
        !is.na(df3_cv$tun_total_CV) & df3_cv$tun_total_CV < cv_cutoff &
        !is.na(df3_cv$thg_total_CV) & df3_cv$thg_total_CV < cv_cutoff,
      na.rm = TRUE
    )
  )
) %>%
  mutate(
    Intensity = factor(Intensity, levels = c("Light", "Heavy", "H/L ratio", "Total"))
  )

p_overlap_after <- ggplot(overlap_after, aes(x = Intensity, y = Protein_Count)) +
  geom_col(width = 0.75, fill = "#6C63B5", color = NA) +
  geom_text(aes(label = Protein_Count), vjust = -0.25, size = 3, color = "black") +
  scale_y_continuous(
    limits = c(0, max(df_facet$Protein_Count, na.rm = TRUE) * 1.12),
    breaks = c(0, 2500, 5000, 7500)
  ) +
  labs(
    x = NULL,
    y = "Protein count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(1, "mm"),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 8)
  )

p_after + p_overlap_after + plot_layout(widths = c(4, 1.2))


#CV boxplot
# 1) reshape CV columns to long
cv_cols <- c(
  "con_light_CV","tun_light_CV","thg_light_CV",
  "con_heavy_CV","tun_heavy_CV","thg_heavy_CV",
  "con_total_CV","tun_total_CV","thg_total_CV",
  "con_HL_CV","tun_HL_CV","thg_HL_CV"
)

cv_long <- df3_cv %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes, all_of(cv_cols)) %>%
  pivot_longer(
    cols = all_of(cv_cols),
    names_to = "key",
    values_to = "CV"
  ) %>%
  filter(is.finite(CV)) %>%
  mutate(
    Treatment = case_when(
      str_detect(key, "^con_") ~ "DMSO",
      str_detect(key, "^tun_") ~ "Tm",
      str_detect(key, "^thg_") ~ "Th",
      TRUE ~ NA_character_
    ),
    Metric = case_when(
      str_detect(key, "_light_") ~ "Light",
      str_detect(key, "_heavy_") ~ "Heavy",
      str_detect(key, "_total_") ~ "Total",
      str_detect(key, "_HL_") ~ "H/L ratio",
      TRUE ~ NA_character_
    ),
    Treatment = factor(Treatment, levels = c("DMSO","Tm","Th")),
    Metric = factor(Metric, levels = c("Light","Heavy","Total","H/L ratio"))
  )

# 2) plot: CV boxplot with 20% dashed line
ggplot(cv_long, aes(x = Treatment, y = CV, fill = Treatment)) +
  geom_boxplot(width = 0.75, color = "black", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", color = "black", shape = 18, size = 2.6) +
  geom_hline(yintercept = 20, linetype = "dashed", linewidth = 0.6) +
  coord_cartesian(ylim = c(0, 50)) +
  scale_fill_manual(values = treat_colors) +
  labs(x = "Treatment", y = "CV (%)") +
  facet_wrap(~Metric, nrow = 1) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.18, "cm"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 11)
  )

#HSPA5 
plot_protein_expression <- function(df, df_name, protein_id) {
  
  
  protein_df <- df %>%
    filter(grepl(protein_id, PG.ProteinAccessions))   
  
  gene_name <- protein_df$PG.Genes[1]
  
  protein_long <- protein_df %>%
    select(
      con_light_mean, tun_light_mean, thg_light_mean,
      con_heavy_mean, tun_heavy_mean, thg_heavy_mean,
      con_total_mean, tun_total_mean, thg_total_mean,
      con_light_SE, tun_light_SE, thg_light_SE,
      con_heavy_SE, tun_heavy_SE, thg_heavy_SE,
      con_total_SE, tun_total_SE, thg_total_SE
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("condition", "channel", "stat"),
      names_pattern = "(con|tun|thg)_(.+?)_(.+)",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = stat,
      values_from = value,
      values_fn = mean
    ) %>%
    mutate(
      channel = factor(channel, 
                       levels = c("light", "heavy", "total"),
                       labels = c("Light", "Heavy", "Total")),
      condition = factor(condition, 
                         levels = c("con","tun", "thg" ),
                         labels = c("DMSO","Tm", "Th"))
    )
  
  
  protein_p_values <- protein_df %>%
    select(light_tun_vs_con_p, light_thg_vs_con_p, heavy_tun_vs_con_p, 
           heavy_thg_vs_con_p, total_tun_vs_con_p, total_thg_vs_con_p)
  
  protein_long <- protein_long %>%
    cross_join(protein_p_values)
  
  
  protein_long <- protein_long %>%
    mutate(
      p_label = case_when(
        condition == "Tm" & channel == "Light" & light_tun_vs_con_p < 0.001 ~ "***",
        condition == "Tm" & channel == "Light" & light_tun_vs_con_p < 0.01 ~ "**",
        condition == "Tm" & channel == "Light" & light_tun_vs_con_p < 0.05 ~ "*",
        
        condition == "Th" & channel == "Light" & light_thg_vs_con_p < 0.001 ~ "***",
        condition == "Th" & channel == "Light" & light_thg_vs_con_p < 0.01 ~ "**",
        condition == "Th" & channel == "Light" & light_thg_vs_con_p < 0.05 ~ "*",
        
        condition == "Tm" & channel == "Heavy" & heavy_tun_vs_con_p < 0.001 ~ "***",
        condition == "Tm" & channel == "Heavy" & heavy_tun_vs_con_p < 0.01 ~ "**",
        condition == "Tm" & channel == "Heavy" & heavy_tun_vs_con_p < 0.05 ~ "*",
        
        condition == "Th" & channel == "Heavy" & heavy_thg_vs_con_p < 0.001 ~ "***",
        condition == "Th" & channel == "Heavy" & heavy_thg_vs_con_p < 0.01 ~ "**",
        condition == "Th" & channel == "Heavy" & heavy_thg_vs_con_p < 0.05 ~ "*",
        
        condition == "Tm" & channel == "Total" & total_tun_vs_con_p < 0.001 ~ "***",
        condition == "Tm" & channel == "Total" & total_tun_vs_con_p < 0.01 ~ "**",
        condition == "Tm" & channel == "Total" & total_tun_vs_con_p < 0.05 ~ "*",
        
        condition == "Th" & channel == "Total" & total_thg_vs_con_p < 0.001 ~ "***",
        condition == "Th" & channel == "Total" & total_thg_vs_con_p < 0.01 ~ "**",
        condition == "Th" & channel == "Total" & total_thg_vs_con_p < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  p <- ggplot(protein_long, aes(x = channel, y = mean, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.75) +
    geom_errorbar(
      aes(ymin = mean - SE, ymax = mean + SE),
      position = position_dodge(width = 0.85),
      width = 0.25
    ) +
    geom_text(
      aes(label = p_label, y = mean + SE + max(mean) * 0.05),
      position = position_dodge(width = 0.85),
      size = 2,
      vjust = 0
    ) +
    scale_y_continuous(
      limits = c(0, 4e6),
      breaks = c(0, 2e6, 4e6)
    ) +
    labs(
      title = paste(gene_name),
      x = "Channel",
      y = "Intensity",
      fill = "Treatment"
    ) +
    scale_fill_manual(values = treat_colors)+
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 8, face = NULL),
      axis.text = element_text(size = 8),
      legend.title =element_blank(),
      legend.text = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black")
    )
  
  return(p)
}
p_df3 <- plot_protein_expression(df3_cv, "df3_cv","P11021")
print(p_df3)

#H/L ratio
plot_HL_expression <- function(df, protein_id) {
  
  protein_df <- df %>%
    filter(grepl(protein_id, PG.ProteinAccessions))
  
  if (nrow(protein_df) == 0) {
    warning("Protein not found.")
    return(NULL)
  }
  
  gene_name <- protein_df$PG.Genes[1]
  
  protein_long <- tibble(
    condition = c("DMSO", "Tm", "Th"),
    mean = c(protein_df$con_HL_mean, protein_df$tun_HL_mean, protein_df$thg_HL_mean),
    SE = c(protein_df$con_HL_SE, protein_df$tun_HL_SE, protein_df$thg_HL_SE),
    p_label = c(
      "",  
      case_when(
        protein_df$HL_tun_vs_con_p < 0.001 ~ "***",
        protein_df$HL_tun_vs_con_p < 0.01 ~ "**",
        protein_df$HL_tun_vs_con_p < 0.05 ~ "*",
        TRUE ~ ""
      ),
      case_when(
        protein_df$HL_thg_vs_con_p < 0.001 ~ "***",
        protein_df$HL_thg_vs_con_p < 0.01 ~ "**",
        protein_df$HL_thg_vs_con_p < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  )
  protein_long$condition <- factor(protein_long$condition, levels = c("DMSO", "Tm", "Th"))
  ggplot(protein_long, aes(x = condition, y = mean, fill = condition)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.2) +
    geom_text(
      aes(label = p_label, y = mean + SE + max(mean, na.rm = TRUE) * 0.05),
      size = 2, vjust = 0
    ) +
    scale_fill_manual(values = treat_colors, name = NULL) +
    scale_y_continuous(
      limits = c(0, max(protein_long$mean + protein_long$SE, na.rm = TRUE) * 1.2),
      breaks = c(0, 2, 4)
    ) +
    labs(
      title = paste0(gene_name),
      y = "log2(H/L ratio)",
      x = NULL
    ) +
    theme_classic() +
    theme(
      plot.title =element_blank(),
      axis.title = element_text(size = 8, face = NULL),
      axis.text = element_text(size = 8),
      legend.title =element_blank(),
      legend.text = element_text(size = 8),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.line = element_line(color = "black")
    )
}
p_df3 <- plot_HL_expression (df3_cv, "P11021")
print(p_df3)


#histogram distribution 
plot_log2_histogram_faceted <- function(df, treatment, binwidth = 0.5) {
  colors <- list("thg" = "#b8b8b8", "tun" = "#ea801c")
  treatment <- tolower(treatment)  
  fill_color <- colors[[treatment]]
  
  
  value_cols <- c("light", "heavy", "HL", "total") %>%
    paste0("_", treatment, "_vs_con")
  
  combined_data <- purrr::map_df(value_cols, function(colname) {
    type <- str_split(colname, "_", simplify = TRUE)[1]
    df_sub <- df %>%
      filter(!is.na(.data[[colname]])) %>%
      mutate(log2_ratio = .data[[colname]],
             Type = case_when(
               type == "light" ~ "Light",
               type == "heavy" ~ "Heavy",
               type == "HL" ~ "H/L",
               type == "total" ~ "Total"
             ))
    return(df_sub)
  })
  combined_data$Type <- factor(combined_data$Type, 
                               levels = c("Light", "Heavy", 
                                          "H/L","Total"))
  annotation_data <- combined_data %>%
    group_by(Type) %>%
    summarise(
      n_count = n(),
      median_val = median(log2_ratio, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(combined_data, aes(x = log2_ratio)) +
    geom_histogram(binwidth = binwidth, fill = fill_color, color = "black") +
    geom_vline(data = annotation_data, aes(xintercept = median_val),
               linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 1) +
    geom_text(data = annotation_data,
              aes(x = 4.5, y = Inf,
                  label = paste0("n = ", n_count, "\nmedian = ", round(median_val, 2))),
              vjust = 1, hjust = 1, size = 2.5) +
    facet_wrap(~Type, nrow = 1, scales = "free_y") +
    scale_x_continuous(limits = c(-4.5, 4.5),breaks = c(-4, -2, 0, 2, 4)) +
    scale_y_continuous(limits = c(0, 5500)) +
    labs(
      title = NULL,
      x = "log2(Ratio)",
      y = "Protein Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      strip.text = element_text(size = 10, color = "black"),
      strip.background = element_rect(fill = "white", color = "black")
    )
}
plot_log2_histogram_faceted(df3_cv, treatment = "Thg")
plot_log2_histogram_faceted(df3_cv, treatment = "Tun")