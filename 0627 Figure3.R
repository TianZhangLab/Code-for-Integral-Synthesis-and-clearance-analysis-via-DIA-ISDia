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

ggplot(df_facet, aes(x = Treatment, y = Protein_Count, fill = Treatment)) +
  geom_col(width = 0.6)+
  geom_text(aes(label = Protein_Count), vjust = -0.3, size = 2, color = "black") +
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
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    theme(panel.spacing = unit(1.5, "cm"))
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

#curve and 48 combine correlation
#all these source df from 0627 Figure 3
df_wide<- df1_Prediction_light %>%
  mutate(log2_AUC = log2(AUC))   %>%
  filter(Condition %in% c("Thg", "DMSO"), Time == 0) %>%
  distinct(PG.ProteinAccessions, PG.Genes, Condition, log2_AUC) %>%
  pivot_wider(names_from = Condition, values_from = log2_AUC) %>%
  drop_na(DMSO, Thg) %>%
  mutate(log2_ratio = Thg- DMSO)
curve_light<-df_wide%>%
  rename(light_thg_vs_con=log2_ratio)


df_wide_heavy<- df1_Prediction %>%
  mutate(log2_AUC = log2(AUC))   %>%
  filter(Condition %in% c("Thg", "DMSO"), Time == 0) %>%
  distinct(PG.ProteinAccessions, PG.Genes, Condition, log2_AUC) %>%
  pivot_wider(names_from = Condition, values_from = log2_AUC) %>%
  drop_na(DMSO, Thg) %>%
  mutate(log2_ratio = Thg- DMSO)
curve_heavy<-df_wide_heavy%>%
  rename(heavy_thg_vs_con=log2_ratio)

df_wide_total<- total_AUC %>%
  mutate(log2_AUC = log2(AUC_total))   %>%
  filter(Condition %in% c("Thg", "DMSO")) %>%
  distinct(PG.ProteinAccessions, PG.Genes, Condition, log2_AUC) %>%
  pivot_wider(names_from = Condition, values_from = log2_AUC) %>%
  drop_na(DMSO, Thg) %>%
  mutate(log2_ratio = Thg- DMSO)
curve_total<-df_wide_total%>%
  rename(total_thg_vs_con=log2_ratio)

df_wide_TR<- df1_Prediction_turnover %>%
  mutate(log2_Turnover_Rate = log2(Turnover_Rate))   %>%
  filter(Treatment %in% c("Thg", "DMSO"), Time == 0) %>%
  distinct(PG.ProteinAccessions, PG.Genes, Treatment, log2_Turnover_Rate) %>%
  pivot_wider(names_from = Treatment, values_from = log2_Turnover_Rate) %>%
  drop_na(DMSO, Thg) %>%
  mutate(log2_ratio = Thg- DMSO)
curve_HL<-df_wide_TR%>%
  rename(HL_thg_vs_con=log2_ratio)
library(purrr)
curve_combined <- list(curve_light, curve_heavy, curve_total,curve_HL) %>%
  purrr::reduce(full_join, by = "PG.ProteinAccessions") %>%
  select(PG.ProteinAccessions, light_thg_vs_con, heavy_thg_vs_con, total_thg_vs_con, HL_thg_vs_con)

combined_WT_48h_curve_combine <- inner_join(curve_combined, df3_cv , by = "PG.ProteinAccessions")

plot_fc_comparison <- function(df, treatment = "thg", channel = "total") {
  fc_col_curve <- paste0(channel, "_", treatment, "_vs_con.x")
  fc_col_combine <- paste0(channel, "_", treatment, "_vs_con.y")
  
  cor_df <- df %>%
    select(PG.ProteinAccessions, 
           curve = all_of(fc_col_curve), 
           combine = all_of(fc_col_combine)) %>%
    filter(!is.na(curve), !is.na(combine))
  
  cor_res <- cor.test(cor_df$curve, cor_df$combine, use = "complete.obs")
  r_val <- round(cor_res$estimate, 2)
  p_val <- formatC(cor_res$p.value, format = "e", digits = 2)
  p_val_display <- ifelse(cor_res$p.value < 2.2e-16, "< 2.2e-16", formatC(cor_res$p.value, format = "e", digits = 2))
  n_val <- df %>%
    filter(!is.na(.data[[fc_col_curve]]), !is.na(.data[[fc_col_combine]])) %>%
    nrow()
  
  p <- ggplot(df, aes(x = .data[[fc_col_curve]], y = .data[[fc_col_combine]])) +
    geom_point(color = "#7ccba2", size = 2, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",  linewidth = 1.2,color = "gray60") +
    scale_x_continuous(limits = c(-5, 5)) +
    scale_y_continuous(limits = c(-5, 5)) +
    annotate("text", x = -4, y = 4,
             label = paste0("R = ", r_val, 
                            "\nP = ", p_val_display, 
                            "\nn = ", n_val),
             size = 5, fontface = "bold", hjust = 0)+
    labs(
      x = paste(str_to_title(channel), treatment, "(48h curve)"),
      y = paste(str_to_title(channel), treatment, "(48h combine)"),
      title = paste(str_to_title(channel), "FC (48h curve vs combine:", treatment, ")")
    )+
    theme_minimal() +
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
  
  return(p)
}
plot_fc_comparison(combined_WT_48h_curve_combine, treatment = "thg", channel = "total")
plot_fc_comparison(combined_WT_48h_curve_combine, treatment = "thg", channel = "light")
plot_fc_comparison(combined_WT_48h_curve_combine, treatment = "thg", channel = "heavy")
plot_fc_comparison(combined_WT_48h_curve_combine, treatment = "thg", channel = "HL")

#The relationship between protein abundance and the regulation of protein synthesis and degradation
df3_separated <- df3_cv %>%
  separate(PG.ProteinAccessions, into = paste0("Accession", 1:25), sep = ";", fill = "right")

PA<-read_excel("D:/personal/UVA/Data analysis/202412 SILAC/abundance_regulation/PA.xlsx")
PA<-PA%>%
  rename(PG.ProteinAccessions="Protein IDs")

PA_separated <- PA %>%
  separate(PG.ProteinAccessions, into = paste0("Accession", 1:45), sep = ";", fill = "right")

df3_merged <- df3_cv %>%
  left_join(PA %>% select(PG.ProteinAccessions, rank), by = "PG.ProteinAccessions")


df3_cv <- df3_cv %>% mutate(RowID = row_number())


df3_long <- df3_cv %>%
  separate(PG.ProteinAccessions, into = paste0("Accession", 1:25), sep = ";", fill = "right") %>%
  pivot_longer(cols = starts_with("Accession"), names_to = "Acc_Index", values_to = "Protein_ID") %>%
  filter(!is.na(Protein_ID))


PA_long <- PA_separated %>%
  pivot_longer(cols = starts_with("Accession"), names_to = "Acc_Index", values_to = "Protein_ID") %>%
  filter(!is.na(Protein_ID)) %>%
  select(Protein_ID, rank) %>%
  distinct()


df3_ranked <- df3_long %>%
  left_join(PA_long, by = "Protein_ID") %>%
  group_by(RowID) %>%
  summarise(rank = first(na.omit(rank)), .groups = "drop")

df3_merged <- df3_cv %>%
  left_join(df3_ranked, by = "RowID") %>%
  select(-RowID)
cat("Number of non-NA ranks:", sum(!is.na(df3_merged$rank)), "\n")


ranked_df <- df3_merged  %>%
  filter(!is.na(rank)) %>%
  arrange(rank) %>%
  mutate(
    rank_group = ntile(rank, 8)  
  )
cat("Number of non-NA ranks:", sum(!is.na(ranked_df$rank)), "\n")
ranked_df%>%
  count(rank_group)


ranked_df_long <- ranked_df %>%
  select(rank_group, light_thg_vs_con, heavy_thg_vs_con, total_thg_vs_con) %>%
  pivot_longer(
    cols = c(light_thg_vs_con, heavy_thg_vs_con, total_thg_vs_con),
    names_to = "Channel",
    values_to = "FC"
  )


rank_colors <- RColorBrewer::brewer.pal(8, "Set2")

#boxplot for LIGHT upregulated proteins
plot_ranked_box_with_stats_up <- function(df, channel,
                                          y_label = NULL,
                                          title = NULL,
                                          colors = rank_colors,
                                          fc_filter_col = "light_thg_vs_con",
                                          pval_filter_col = "light_thg_vs_con_p",
                                          fc_threshold = 0.58,
                                          pval_threshold = 0.05) {
  channel_sym <- sym(channel)
  fc_sym <- sym(fc_filter_col)
  pval_sym <- sym(pval_filter_col)
  
  df_filtered <- df %>%
    filter(
      !is.na(rank_group),
      !!fc_sym > fc_threshold,
      !!pval_sym < pval_threshold
    )
  
  if (nrow(df_filtered) == 0) {
    message("No data points meet the filter criteria.")
    return(NULL)
  }
  
  # stats per group
  stat_df <- df_filtered %>%
    group_by(rank_group = as.factor(rank_group)) %>%
    summarise(
      mean_val = mean(!!channel_sym, na.rm = TRUE),
      median_val = median(!!channel_sym, na.rm = TRUE),
      sd_val = sd(!!channel_sym, na.rm = TRUE),
      n = sum(!is.na(!!channel_sym)),
      .groups = "drop"
    )
  
  # t-tests against group 1 (one-tailed, unpaired)
  group1_vals <- df_filtered %>%
    filter(rank_group == 1) %>%
    pull(!!channel_sym)
  
  pval_df <- df_filtered %>%
    filter(rank_group != 1) %>%
    group_by(rank_group = as.factor(rank_group)) %>%
    summarise(
      p_value = tryCatch({
        group_vals <- pull(cur_data(), !!channel_sym)
        t.test(group1_vals, group_vals, alternative = "less")$p.value
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value > 0.05  ~ "ns",
      TRUE            ~ ""
    ))
  
  stat_df <- stat_df %>%
    left_join(pval_df, by = "rank_group")
  
  p <- ggplot(df_filtered, aes(x = as.factor(rank_group), y = !!channel_sym, fill = as.factor(rank_group))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_text(data = stat_df,
              aes(
                x = rank_group,
                y = 0.25,
                label = paste0("n=", n, "\nMed=", round(median_val, 2), "\nSD=", round(sd_val, 2))
              ),
              size = 2.5, color = "black", inherit.aes = FALSE) +
    geom_text(data = stat_df,
              aes(
                x = rank_group,
                y = 1.7 ,
                label = stars
              ),
              size = 4, color = "red", inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0,1.8))+
    scale_fill_manual(values = colors) +
    labs(
      title = title %||% paste(channel, "by Rank Group"),
      x = "Rank Group",
      y = y_label %||% channel
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "none" )
  
  return(p)
}
plot_ranked_box_with_stats_up(ranked_df, "light_thg_vs_con", title = "Light FC by Rank Group")


#boxplot for heavy upregulated proteins
plot_ranked_box_with_stats_up <- function(df, channel,
                                          y_label = NULL,
                                          title = NULL,
                                          colors = rank_colors,
                                          fc_filter_col = "heavy_thg_vs_con",
                                          pval_filter_col = "heavy_thg_vs_con_p",
                                          fc_threshold = 0.58,
                                          pval_threshold = 0.05) {
  channel_sym <- sym(channel)
  fc_sym <- sym(fc_filter_col)
  pval_sym <- sym(pval_filter_col)
  
  df_filtered <- df %>%
    filter(
      !is.na(rank_group),
      !!fc_sym > fc_threshold,
      !!pval_sym < pval_threshold
    )
  
  if (nrow(df_filtered) == 0) {
    message("No data points meet the filter criteria.")
    return(NULL)
  }
  
  # stats per group
  stat_df <- df_filtered %>%
    group_by(rank_group = as.factor(rank_group)) %>%
    summarise(
      mean_val = mean(!!channel_sym, na.rm = TRUE),
      median_val = median(!!channel_sym, na.rm = TRUE),
      sd_val = sd(!!channel_sym, na.rm = TRUE),
      n = sum(!is.na(!!channel_sym)),
      .groups = "drop"
    )
  
  # t-tests against group 1 (one-tailed, unpaired)
  group1_vals <- df_filtered %>%
    filter(rank_group == 1) %>%
    pull(!!channel_sym)
  
  pval_df <- df_filtered %>%
    filter(rank_group != 1) %>%
    group_by(rank_group = as.factor(rank_group)) %>%
    summarise(
      p_value = tryCatch({
        group_vals <- pull(cur_data(), !!channel_sym)
        t.test(group1_vals, group_vals, alternative = "less")$p.value
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value > 0.05  ~ "ns",
      TRUE            ~ ""
    ))
  
  stat_df <- stat_df %>%
    left_join(pval_df, by = "rank_group")

  p <- ggplot(df_filtered, aes(x = as.factor(rank_group), y = !!channel_sym, fill = as.factor(rank_group))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_text(data = stat_df,
              aes(
                x = rank_group,
                y = 0.25,
                label = paste0("n=", n, "\nMed=", round(median_val, 2), "\nSD=", round(sd_val, 2))
              ),
              size = 2.5, color = "black", inherit.aes = FALSE) +
    geom_text(data = stat_df,
              aes(
                x = rank_group,
                y = 2.5 ,
                label = stars
              ),
              size = 4, color = "red", inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0,2.6))+
    scale_fill_manual(values = colors) +
    labs(
      title = title %||% paste(channel, "by Rank Group"),
      x = "Rank Group",
      y = y_label %||% channel
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "none" ) 
  
  return(p)
}
plot_ranked_box_with_stats_up(ranked_df, "heavy_thg_vs_con", title = "heavy FC by Rank Group")


