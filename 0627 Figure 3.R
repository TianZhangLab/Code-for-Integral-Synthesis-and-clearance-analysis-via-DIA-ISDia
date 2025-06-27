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
library(VennDiagram)

#df1=curve,df2=48h alone, df3=combine timepoint
df1 <- read.delim("D:/personal/UVA/Data analysis/202412 SILAC/A012_treated/20250204_132507_A012_SILAC_ALL2_Report with gene.tsv", header = TRUE, sep = "\t")

#Curve method
#df1 CV calculate and count
transform_data <- function(df) {
  ms2_cols <- grep("PG\\.MS2Channel[12]", colnames(df), value = TRUE)
  df_long <- df %>%
    pivot_longer(cols = all_of(ms2_cols), names_to = "Sample", values_to = "Intensity") %>%
    mutate(
      # Assign Time based on Sample name
      Time = case_when(
        str_detect(Sample, "thg_[123]\\.") | str_detect(Sample, "tun_[123]\\.") | str_detect(Sample, "con_[456]\\.") ~ "6h",
        str_detect(Sample, "thg_[456]\\.") | str_detect(Sample, "tun_[456]\\.") | str_detect(Sample, "con_[789]\\.") ~ "12h",
        str_detect(Sample, "thg_[789]\\.") | str_detect(Sample, "tun_[789]\\.") | str_detect(Sample, "con_1[012]\\.") ~ "24h",
        str_detect(Sample, "thg_1[012]\\.") | str_detect(Sample, "tun_1[012]\\.") | str_detect(Sample, "con_1[345]\\.") ~ "48h",
        TRUE ~ NA_character_
      ),
      # Assign Type based on MS2Channel
      Type = case_when(
        str_detect(Sample, "PG\\.MS2Channel1") ~ "Light",
        str_detect(Sample, "PG\\.MS2Channel2") ~ "Heavy",
        TRUE ~ NA_character_
      ),
      Treatment=case_when(
        str_detect(Sample, "thg_\\d.") ~ "Thg",
        str_detect(Sample,"tun_\\d.") ~ "Tun",
        str_detect(Sample,"con_\\d.") ~ "DMSO",
        TRUE ~ NA_character_
      )
    )
  # Check if any Time or Type is NA and throw an error
  if (any(is.na(df_long$Time))) {
    warning("Some samples did not match the expected patterns ('thg_', 'tun_', or 'con_'). Check your column names.")
  }
  if (any(is.na(df_long$Type))) {
    warning("Some samples did not match 'PG.MS2Channel1' or 'PG.MS2Channel2'. Check your column names.")
  }
  
  return(df_long)
}

#apply this function to df
df1_transformed <- transform_data(df1)

# calculate_cv
calculate_cv <- function(df_long) {
  df_long %>%
    group_by(PG.ProteinAccessions, Treatment, Time, Type) %>%
    summarise(
      CV = ifelse(
        sum(!is.na(as.numeric(Intensity))) >= 3,  
        sd(as.numeric(Intensity), na.rm = TRUE) / mean(as.numeric(Intensity), na.rm = TRUE) * 100,  
        NA_real_  
      ),
      Mean = ifelse(
        sum(!is.na(as.numeric(Intensity))) >= 3,
        mean(as.numeric(Intensity), na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = c(Time,Treatment, Type), values_from =c(CV, Mean),names_sep = "_")
}
df1_cv <- calculate_cv(df1_transformed)
df1_cv <- df1 %>%
  left_join(df1_cv, by = "PG.ProteinAccessions")

get_time_label <- function(sample_name) {
  case_when(
    # 6h conditions
    str_detect(sample_name, "thg_[123]\\.") |
      str_detect(sample_name, "tun_[123]\\.") |
      str_detect(sample_name, "con_[456]\\.") ~ "6h",
    
    # 12h conditions
    str_detect(sample_name, "thg_[456]\\.") |
      str_detect(sample_name, "tun_[456]\\.") |
      str_detect(sample_name, "con_[789]\\.") ~ "12h",
    
    # 24h conditions
    str_detect(sample_name, "thg_[789]\\.") |
      str_detect(sample_name, "tun_[789]\\.") |
      str_detect(sample_name, "con_1[012]\\.") ~ "24h",
    
    # 48h conditions
    str_detect(sample_name, "thg_1[012]\\.") |
      str_detect(sample_name, "tun_1[012]\\.") |
      str_detect(sample_name, "con_1[345]\\.") ~ "48h",
    
    # Default
    TRUE ~ NA_character_
  )
}

#2 time point synthesis curve with real data
# ---- Start of code block ----
clean_and_transform_heavy <- function(df) {
  df_heavy <- df %>%
    select(-contains("Light"))
  
  intensity_cols <- colnames(df_heavy)[str_detect(colnames(df_heavy), "MS2Channel2")]
  for (col in intensity_cols) {
    time <- get_time_label(col)
    treat <- case_when(
      str_detect(col, "_con_") ~ "DMSO",
      str_detect(col, "_tun_") ~ "Tun",
      str_detect(col, "_thg_") ~ "Thg",
      TRUE ~ NA_character_
    )
    cv_col <- paste0("CV_", time, "_", treat, "_Heavy")
    if (!is.na(time) && !is.na(treat) && cv_col %in% colnames(df_heavy)) {
      df_heavy[[col]] <- ifelse(df_heavy[[cv_col]] > 20, NA, df_heavy[[col]])
    }
  }
  
  df_long_heavy <- df_heavy %>%
    select(PG.ProteinAccessions, PG.Genes,all_of(intensity_cols)) %>%
    pivot_longer(cols = all_of(intensity_cols),
                 names_to = "raw_colname",
                 values_to = "intensity") %>%
    filter(!is.na(intensity)) %>%
    mutate(
      treatment = case_when(
        str_detect(raw_colname, "_con_") ~ "DMSO",
        str_detect(raw_colname, "_tun_") ~ "Tun",
        str_detect(raw_colname, "_thg_") ~ "Thg",
        TRUE ~ NA_character_
      ),
      time = get_time_label(raw_colname),
      Time = as.numeric(str_remove(time, "h")),
      Condition = treatment
    ) %>%
    filter(!is.na(treatment), !is.na(time)) %>%
    select(PG.ProteinAccessions, PG.Genes,Intensity = intensity, Time, Condition)
  
  df_long_heavy_filtered <- df_long_heavy %>%
    group_by(PG.ProteinAccessions, Condition) %>%
    filter(n() >= 6) %>%
    ungroup()
  
  return(df_long_heavy_filtered)
}
filtered_df1 <- clean_and_transform_heavy(df1_cv)
length(unique(filtered_df1$PG.ProteinAccessions))

synthesis_curve_2TP_real <- function(df_with_cv) {
  df_with_cv <- bind_rows(df_with_cv,
                          df_with_cv %>%
                            group_by(PG.ProteinAccessions, Condition) %>%
                            slice(1) %>%
                            ungroup() %>%
                            mutate(Time = 0, Intensity = 0)
  ) 
  
  synthesis_model <- function(time, S, k) {
    S * (1 - exp(-k * time))
  }
  
  fit_synthesis_model <- function(df) {
    tryCatch({
      fit <- nlsLM(Intensity ~ synthesis_model(Time, S, k),
                   data = df,
                   start = list(S = max(df$Intensity, na.rm = TRUE), k = 0.1),
                   control = nls.lm.control(maxiter = 500))
      
      coef_vals <- coef(summary(fit))
      S <- coef_vals["S", "Estimate"]
      k <- coef_vals["k", "Estimate"]
      S_SE <- coef_vals["S", "Std. Error"]
      k_SE <- coef_vals["k", "Std. Error"]
      k_t <- coef_vals["k", "t value"]
      k_p <- coef_vals["k", "Pr(>|t|)"]
      S_p <- coef_vals["S", "Pr(>|t|)"]
      auc <- S * (50 + (1 / k) * (exp(-k * 50) - 1))
      pred <- predict(fit)
      residuals <- residuals(fit)
      R_squared <- 1 - (sum(residuals^2) / sum((df$Intensity - mean(df$Intensity))^2))
      Residual_SD <- sqrt(sum(residuals(fit)^2) / df.residual(fit))
      Residual_SD_S_Ratio <- Residual_SD / S
      
      return(data.frame(
        PG.ProteinAccessions = unique(df$PG.ProteinAccessions),
        Condition = unique(df$Condition),
        Synthesis_constant = S,
        Synthesis_Rate = k,
        S_SE = S_SE,
        k_SE = k_SE,
        k_tvalue = k_t,
        k_pvalue = k_p,
        S_pvalue = S_p,
        Residual_SD = Residual_SD,
        R_squared = R_squared,
        Residual_SD_S_Ratio=Residual_SD_S_Ratio,
        AUC = auc
      ))
    }, error = function(e) {
      return(data.frame(
        PG.ProteinAccessions = unique(df$PG.ProteinAccessions),
        Condition = unique(df$Condition),
        Synthesis_constant = NA,
        Synthesis_Rate = NA,
        S_SE = NA,
        k_SE = NA,
        k_tvalue = NA,
        k_pvalue = NA,
        S_pvalue = NA,
        Residual_SD = NA,
        R_squared = NA,
        Residual_SD_S_Ratio=NA,
        AUC = NA_real_
      ))
    })
  }
  
  protein_k_values <- df_with_cv %>%
    group_by(PG.ProteinAccessions, Condition) %>%
    group_split() %>%
    map_df(fit_synthesis_model)
  
  df_with_cv <- df_with_cv %>%
    left_join(protein_k_values, by = c("PG.ProteinAccessions", "Condition"))
  
  time_seq <- seq(0, 50, by = 1)
  prediction_data <- expand.grid(
    PG.ProteinAccessions = unique(protein_k_values$PG.ProteinAccessions),
    Condition = unique(protein_k_values$Condition),  
    Time = time_seq
  ) %>%
    left_join(protein_k_values, by = c("PG.ProteinAccessions", "Condition")) %>%
    mutate(Predicted_Intensity = synthesis_model(Time, Synthesis_constant, Synthesis_Rate)) %>%
    left_join(dplyr::select(df_with_cv, PG.ProteinAccessions, PG.Genes, Condition) %>% distinct(), 
              by = c("PG.ProteinAccessions", "Condition")) 
  
  return(list(
    Long_Filtered_Data = df_with_cv,
    Prediction_Data = prediction_data,
    R_squared = protein_k_values
  ))
}
result_df1_2TP_real <- synthesis_curve_2TP_real(filtered_df1)
df1_Prediction <- result_df1_2TP_real$Prediction_Data
df1_filtered <- result_df1_2TP_real$Long_Filtered_Data
length(unique(df1_Prediction$PG.ProteinAccessions))


#SINGLE PROTEIN synthesis CURVE
plot_synthesis_curve <- function(df, df_name, protein_id, model_info_df = NULL, raw_points = NULL) {
  df_filtered <- df %>% 
    filter(PG.ProteinAccessions == protein_id)%>%
    mutate(Condition = recode(Condition,
                              "Tun" = "Tm",
                              "Thg" = "Th"))

  condition_colors <- c("DMSO" = "#1a80bb", "Tm" = "#ea801c", "Th" = "#b8b8b8")
  
  p <- ggplot(df_filtered, aes(x = Time, y = Predicted_Intensity, color = Condition)) +
    scale_color_manual(values = condition_colors) +
    scale_y_continuous(labels = scales::scientific)+
    geom_line(linewidth = 1) +
    labs(title = df_filtered$PG.Genes[which(!is.na(df_filtered$PG.Genes))[1]],
         x = "Time (hours)",
         y = "Intensity") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = NULL),
      axis.text = element_text(size = 12),
      legend.title =element_blank(),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "bottom"
    )
  
 
  if (!is.null(raw_points)) {
    raw_subset <- raw_points %>% filter(PG.ProteinAccessions == protein_id)%>%
      mutate(Condition = recode(Condition,
                                "Tun" = "Tm",
                                "Thg" = "Th"))
    p <- p + geom_point(data = raw_subset, aes(x = Time, y = Intensity, color = Condition), shape = 16, size = 2)
  }
  
 
  if (!is.null(model_info_df)) {
    model_subset <- model_info_df %>%
      filter(PG.ProteinAccessions == protein_id) %>%
      mutate(Condition = recode(Condition, "Tun" = "Tm", "Thg" = "Th"))  # ⚠️ 也改名
    
    if (nrow(model_subset) > 0) {
      label_text <- model_subset %>%
        mutate(label = paste0(Condition, ": R²=", round(R_squared, 2))) %>%
        pull(label) %>%
        paste(collapse = "\n")
      
      p <- p + annotate("text",
                        x = 5,
                        y = max(df_filtered$Predicted_Intensity, na.rm = TRUE) * 0.8,
                        label = label_text,
                        hjust = 0,
                        size = 2.5,
                        color = "black")
    }
  }
  
  return(p)
}
p_df1_heavy <- plot_synthesis_curve(df = df1_Prediction,df_name = "df1",protein_id = "P11021",model_info_df = result_df1_2TP_real$R_squared,raw_points = filtered_df1) 
print(p_df1_heavy)
# ---- End of code block ----



#2 time point degradation curve with real data
# ---- Start of code block ----
clean_and_transform_light <- function(df) {
  df_light <- df %>%
    select(-contains("Heavy"))
  
  intensity_cols <- colnames(df_light)[str_detect(colnames(df_light), "MS2Channel1")]
  
  for (col in intensity_cols) {
    time <- get_time_label(col)
    treat <- case_when(
      str_detect(col, "_con_") ~ "DMSO",
      str_detect(col, "_tun_") ~ "Tun",
      str_detect(col, "_thg_") ~ "Thg",
      TRUE ~ NA_character_
    )
    
    cv_col <- paste0("CV_", time, "_", treat, "_Light")
    
    if (!is.na(time) && !is.na(treat) && cv_col %in% colnames(df_light)) {
      df_light[[col]] <- ifelse(df_light[[cv_col]] > 20, NA, df_light[[col]])
    }
  }
  
  df_long <- df_light %>%
    select(PG.ProteinAccessions, PG.Genes,all_of(intensity_cols)) %>%
    pivot_longer(cols = all_of(intensity_cols),
                 names_to = "raw_colname",
                 values_to = "intensity") %>%
    filter(!is.na(intensity)) %>%
    mutate(
      treatment = case_when(
        str_detect(raw_colname, "_con_") ~ "DMSO",
        str_detect(raw_colname, "_tun_") ~ "Tun",
        str_detect(raw_colname, "_thg_") ~ "Thg",
        TRUE ~ NA_character_
      ),
      time = get_time_label(raw_colname)
    ) %>%
    filter(!is.na(treatment), !is.na(time))
  
  df_long_filtered <- df_long %>%
    group_by(PG.ProteinAccessions, treatment) %>%
    filter(n() >= 6) %>%
    ungroup()
  
  final_df <- df_long_filtered %>%
    mutate(
      Time = as.numeric(str_remove(time, "h")),
      Condition = treatment
    ) %>%
    select(PG.ProteinAccessions, PG.Genes,Intensity = intensity, Time, Condition)
  
  return(final_df)
}
filtered_df1_light <- clean_and_transform_light(df1_cv)
length(unique(filtered_df1_light$PG.ProteinAccessions))

#set I0 as the highest intensity for the same treatment in same protein,but for same protein with different conditions, they have individual I0
degradation_curve_2TP_real<- function(df_with_cv) {
  
  fit_individual_I0 <- function(df) {
    tryCatch({
      fit <- nlsLM(Intensity ~ I0 * exp(-k * Time),
                   data = df,
                   start = list(I0 = max(df$Intensity, na.rm = TRUE), k = 0.1),
                   control = nls.lm.control(maxiter = 500))
      I0 <- coef(fit)[["I0"]]
      
      if (any(df$Intensity >= I0, na.rm = TRUE)) {
        I0 <- max(df$Intensity, na.rm = TRUE) * 1.01  
      }
      tibble(PG.ProteinAccessions = unique(df$PG.ProteinAccessions), 
             Condition = unique(df$Condition),
             I0 = I0)
    }, error = function(e) {
      tibble(PG.ProteinAccessions = unique(df$PG.ProteinAccessions), 
             Condition = unique(df$Condition),
             I0 = NA_real_)
    })
  }
  
  i0_per_protein_condition <- df_with_cv %>%
    group_by(PG.ProteinAccessions, Condition) %>%
    group_split() %>%
    map_df(fit_individual_I0)
  
  long_with_i0 <- df_with_cv %>%
    left_join(i0_per_protein_condition, by = c("PG.ProteinAccessions", "Condition")) %>%
    filter(!is.na(I0))
  
  fit_k_fixed_I0 <- function(df) {
    tryCatch({
      I0_value <- unique(df$I0)
      fit <- nlsLM(Intensity ~ I0_value * exp(-k * Time),
                   data = df,
                   start = list(k = 0.1),
                   control = nls.lm.control(maxiter = 500))
      k <- coef(fit)[["k"]]
      coef_vals <- coef(summary(fit)) 
      k_SE <- coef_vals["k", "Std. Error"]
      k_t <- coef_vals["k", "t value"]
      k_p <- coef_vals["k", "Pr(>|t|)"]
      residuals <- residuals(fit)
      auc <- (I0_value / k) * (1 - exp(-k * 50))
      R_squared <- 1 - (sum(residuals^2) / sum((df$Intensity - mean(df$Intensity))^2))
      Residual_SD <- sqrt(sum(residuals^2) / df.residual(fit))
      Residual_SD_I0_Ratio <- Residual_SD / I0_value
      
      tibble(
        PG.ProteinAccessions = unique(df$PG.ProteinAccessions),
        Condition = unique(df$Condition),
        Initial_Intensity = I0_value,
        Degradation_Rate = k,
        R_squared = R_squared,
        Residual_SD = Residual_SD,
        Residual_SD_I0_Ratio = Residual_SD_I0_Ratio,
        AUC = auc,
        k_SE = k_SE,
        k_tvalue = k_t,
        k_pvalue = k_p
      )
    }, error = function(e) {
      tibble(
        PG.ProteinAccessions = unique(df$PG.ProteinAccessions),
        Condition = unique(df$Condition),
        Initial_Intensity = NA,
        Degradation_Rate = NA,
        R_squared = NA,
        Residual_SD = NA,
        Residual_SD_I0_Ratio = NA,
        AUC = NA,
        k_SE = NA,
        k_tvalue = NA,
        k_pvalue = NA
      )
    })
  }
  
  protein_degradation_rates <- long_with_i0 %>%
    group_by(PG.ProteinAccessions, Condition) %>%
    group_split() %>%
    map_df(fit_k_fixed_I0)
  
  degradation_model <- function(time, I0, k) {
    I0 * exp(-k * time)
  }
  
  time_seq <- seq(0, 50, by = 1)
  prediction_data <- expand.grid(
    PG.ProteinAccessions = unique(protein_degradation_rates$PG.ProteinAccessions),
    Condition = unique(protein_degradation_rates$Condition),
    Time = time_seq
  ) %>%
    left_join(protein_degradation_rates, by = c("PG.ProteinAccessions", "Condition")) %>%
    mutate(Predicted_Intensity = degradation_model(Time, Initial_Intensity, Degradation_Rate))%>%
    left_join(dplyr::select(df_with_cv, PG.ProteinAccessions, PG.Genes, Condition) %>% distinct(), 
              by = c("PG.ProteinAccessions", "Condition")) 
  
  return(list(
    Protein_Degradation_Rates = protein_degradation_rates,
    Prediction_Data = prediction_data
  ))
}
result_df1_2TP_real_light <- degradation_curve_2TP_real(filtered_df1_light)
df1_Prediction_light <- result_df1_2TP_real_light$Prediction_Data
length(unique(df1_Prediction_light$PG.ProteinAccessions))

#SINGLE PROTEIN DEGRADATION CURVE
plot_degradation_curve <- function(df, df_name, protein_id, model_info_df = NULL, raw_points = NULL) {
  df_filtered <- df %>% filter(PG.ProteinAccessions == protein_id)%>%
    mutate(Condition = recode(Condition,
                              "Tun" = "Tm",
                              "Thg" = "Th"))

  condition_colors <- c("DMSO" = "#1a80bb", "Tm" = "#ea801c", "Th" = "#b8b8b8")
  
  p <- ggplot(df_filtered, aes(x = Time, y = Predicted_Intensity, color = Condition)) +
    scale_color_manual(values = condition_colors) +
    scale_y_continuous(labels = scales::scientific)+
    geom_line(linewidth = 1) +
    labs(title = df_filtered$PG.Genes[which(!is.na(df_filtered$PG.Genes))[1]],
         x = "Time (hours)",
         y = "Intensity") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = NULL),
      axis.text = element_text(size = 12),
      legend.title =element_blank(),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "bottom"
    )
  
  if (!is.null(raw_points)) {
    raw_subset <- raw_points %>% filter(PG.ProteinAccessions == protein_id)%>%
      mutate(Condition = recode(Condition,
                                "Tun" = "Tm",
                                "Thg" = "Th"))
    p <- p + geom_point(data = raw_subset, aes(x = Time, y = Intensity, color = Condition), shape = 16, size = 2.5)
  }
  

  if (!is.null(model_info_df)) {
    model_subset <- model_info_df %>%
      filter(PG.ProteinAccessions == protein_id) %>%
      mutate(Condition = recode(Condition, "Tun" = "Tm", "Thg" = "Th"))  
    
    if (nrow(model_subset) > 0) {
      label_text <- model_subset %>%
        mutate(label = paste0(Condition, ": R²=", round(R_squared, 2))) %>%
        pull(label) %>%
        paste(collapse = "\n")
      
      p <- p + annotate("text",
                        x = 35,
                        y = max(df_filtered$Predicted_Intensity, na.rm = TRUE) * 0.8,
                        label = label_text,
                        hjust = 0,
                        size = 2.5,
                        color = "black")
    }
  }
  return(p)
}
p_df1_light <- plot_degradation_curve(df = df1_Prediction_light,df_name = "df1",protein_id = "P11021",model_info_df = result_df1_2TP_real_light$Protein_Degradation_Rates,raw_points = filtered_df1_light) 
print(p_df1_light)
# ---- End of code block ----

#Total AUC by predict sum of light and heavy AUC
heavy_auc <- df1_Prediction %>% filter(Time == 48) %>%
  select(PG.ProteinAccessions, Condition, AUC_heavy = AUC,PG.Genes)
light_auc <- df1_Prediction_light %>% filter(Time == 0) %>%
  select(PG.ProteinAccessions, Condition, AUC_light = AUC)
merged_df <- inner_join(heavy_auc, light_auc, by = c("PG.ProteinAccessions", "Condition"))
total_AUC <- inner_join(heavy_auc, light_auc, by = c("PG.ProteinAccessions", "Condition")) %>%
  filter(!is.na(AUC_heavy), !is.na(AUC_light)) %>%
  mutate(AUC_total = AUC_heavy + AUC_light)

#single protein AUC comparision plot
plot_protein_auc_by_type <- function(df, protein_id) {
  
  treat_colors <- c("DMSO" = "#1a80bb", "Tm" = "#ea801c", "Th" = "#b8b8b8")
  
  protein_df <- df %>%
    filter(str_detect(PG.ProteinAccessions, protein_id)) 
  
  if (nrow(protein_df) == 0) {
    stop("No matching protein found.")
  }
  
  gene_name <- unique(protein_df$PG.Genes)[1]
  
  protein_long <- protein_df %>%
    pivot_longer(cols = c(AUC_light, AUC_heavy, AUC_total),
                 names_to = "Type", values_to = "AUC") %>%
    mutate(
      Type = recode(Type,
                    AUC_light = "Light",
                    AUC_heavy = "Heavy",
                    AUC_total = "Total"),
      Type = factor(Type, levels = c("Light", "Heavy", "Total")),
      Condition = recode(Condition, "Tun" = "Tm", "Thg" = "Th"),
      Condition = factor(Condition, levels = c("DMSO", "Tm", "Th"))
    )
  

  ggplot(protein_long, aes(x = Type, y = AUC, fill = Condition)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.6, color = "black") +
    labs(
      title = gene_name,
      x = NULL,
      y = "AUC",
      fill = "Condition"
    ) +
    scale_fill_manual(values = treat_colors) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 14)
    )
}
plot_protein_auc_by_type(total_AUC , protein_id = "P11021")


#Turnover calculation and curve fitting
# ---- Start of code block ----
df1_transformed<- df1_transformed %>%
  mutate(
    replicate = case_when(
      str_detect(Sample, "thg_(1[0]|[147])\\.raw") |
        str_detect(Sample, "tun_(1[0]|[147])\\.raw") |
        str_detect(Sample, "con_(1[03]|[47])\\.raw") ~ "1",
      
      str_detect(Sample, "thg_(1[1]|[258])\\.raw") |
        str_detect(Sample, "tun_(1[1]|[258])\\.raw") |
        str_detect(Sample, "con_(1[14]|[58])\\.raw") ~ "2",
      
      str_detect(Sample, "thg_(1[2]|[369])\\.raw") |
        str_detect(Sample, "tun_(1[2]|[369])\\.raw") |
        str_detect(Sample, "con_(1[25]|[69])\\.raw") ~ "3",
      
      TRUE ~ NA_character_
    )
  )

LT_ratio<- df1_transformed %>%
  filter(!is.na(Intensity)) %>%
  pivot_wider(
    id_cols = c(PG.ProteinAccessions, PG.Genes, Treatment, Time,replicate),
    names_from = Type,
    values_from = Intensity
  ) %>%
  mutate(
    LT_ratio = ifelse(!is.na(Light) & !is.na(Heavy), Light / (Light + Heavy), NA_real_)
  ) %>%
  group_by(PG.ProteinAccessions, Treatment, Time) %>%
  mutate(
    LT_ratio_CV = {
      values <- LT_ratio[!is.na(LT_ratio)]
      if (length(values) >= 3) {
        (sd(values) / mean(values)) * 100
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()


#filter protein with at least one time point CV<20% for fitting
LT_ratio_filtered <- LT_ratio %>%
  filter(!is.na(LT_ratio_CV), LT_ratio_CV <= 20) %>%
  group_by(PG.ProteinAccessions, Treatment) %>%
  filter(n_distinct(Time) >= 1) %>%
  mutate(Time = as.numeric(stringr::str_remove(Time, "h")))%>%
  ungroup()


#fit curve
turnover_curve <- function(df_with_cv) {
  
  # Add(0,1) 
  df_with_cv <- bind_rows(
    df_with_cv,
    df_with_cv %>%
      group_by(PG.ProteinAccessions, Treatment) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(Time = 0, LT_ratio = 1)
  )
  
  # LT_ratio = exp(-k * time)
  turnover_model <- function(time, k) {
    exp(-k * time)
  }
  
  fit_turnover_model <- function(df) {
    tryCatch({
      fit <- nlsLM(
        LT_ratio ~ turnover_model(Time, k),
        data = df,
        start = list(k = 0.1),
        control = nls.lm.control(maxiter = 500)
      )
      
      coef_vals <- coef(summary(fit))
      k <- coef_vals["k", "Estimate"]
      k_SE <- coef_vals["k", "Std. Error"]
      k_t <- coef_vals["k", "t value"]
      k_p <- coef_vals["k", "Pr(>|t|)"]
      
      residuals <- residuals(fit)
      R_squared <- 1 - (sum(residuals^2) / sum((df$LT_ratio - mean(df$LT_ratio))^2))
      Residual_SD <- sqrt(sum(residuals^2) / df.residual(fit))
      
      tibble(
        PG.ProteinAccessions = unique(df$PG.ProteinAccessions),
        Treatment = unique(df$Treatment),
        Turnover_Rate = k,
        k_SE = k_SE,
        k_tvalue = k_t,
        k_pvalue = k_p,
        Residual_SD = Residual_SD,
        R_squared = R_squared
      )
    }, error = function(e) {
      tibble(
        PG.ProteinAccessions = unique(df$PG.ProteinAccessions),
        Treatment = unique(df$Treatment),
        Turnover_Rate = NA_real_,
        k_SE = NA_real_,
        k_tvalue = NA_real_,
        k_pvalue = NA_real_,
        Residual_SD = NA_real_,
        R_squared = NA_real_
      )
    })
  }
  
  protein_k_values <- df_with_cv %>%
    group_by(PG.ProteinAccessions, Treatment) %>%
    group_split() %>%
    map_df(fit_turnover_model)
  
  df_with_cv <- df_with_cv %>%
    left_join(protein_k_values, by = c("PG.ProteinAccessions", "Treatment"))
  
  time_seq <- seq(0, 50, by = 1)
  prediction_data <- expand.grid(
    PG.ProteinAccessions = unique(protein_k_values$PG.ProteinAccessions),
    Treatment = unique(protein_k_values$Treatment),
    Time = time_seq
  ) %>%
    left_join(protein_k_values, by = c("PG.ProteinAccessions", "Treatment")) %>%
    mutate(Predicted_LT_ratio = turnover_model(Time, Turnover_Rate)) %>%
    left_join(
      df_with_cv %>%
        select(PG.ProteinAccessions, PG.Genes, Treatment) %>%
        distinct(),
      by = c("PG.ProteinAccessions", "Treatment")
    )
  
  return(list(
    LT_ratio_with_fits = df_with_cv,
    Prediction_Data = prediction_data,
    Protein_turnover_Rates = protein_k_values
  ))
}
result_df1_turnover<- turnover_curve (LT_ratio_filtered)
df1_Prediction_turnover <- result_df1_turnover$Prediction_Data

#SINGLE PROTEIN turnover CURVE
plot_turnover_curve <- function(df, df_name, protein_id, model_info_df = NULL, raw_points = NULL) {
  
  df_filtered <- df %>% 
    filter(PG.ProteinAccessions == protein_id) %>%
    mutate(Condition = recode(Treatment,
                              "Tun" = "Tm",
                              "Thg" = "Th"))
  
  condition_colors <- c("DMSO" = "#1a80bb", "Tm" = "#ea801c", "Th" = "#b8b8b8")
  p <- ggplot(df_filtered, aes(x = Time, y = Predicted_LT_ratio, color = Condition)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = condition_colors) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = df_filtered$PG.Genes[which(!is.na(df_filtered$PG.Genes))[1]],
      x = "Time (hours)",
      y = "Light/Total ratio"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "bottom"
    )
  
  if (!is.null(raw_points)) {
    raw_subset <- raw_points %>%
      filter(PG.ProteinAccessions == protein_id) %>%
      mutate(Condition = recode(Treatment,
                                "Tun" = "Tm",
                                "Thg" = "Th"))
    
    p <- p + geom_point(data = raw_subset,
                        aes(x = Time, y = LT_ratio, color = Condition),
                        shape = 16, size = 2.5)
  }
  
  if (!is.null(model_info_df)) {
    model_subset <- model_info_df %>%
      filter(PG.ProteinAccessions == protein_id) %>%
      mutate(Condition = recode(Treatment, "Tun" = "Tm", "Thg" = "Th"))
    
    if (nrow(model_subset) > 0) {
      label_text <- model_subset %>%
        mutate(
          label = paste0(
            Condition,
            ": k=", round(Turnover_Rate, 3),
            ", R²=", round(R_squared, 2)
          )
        ) %>%
        pull(label) %>%
        paste(collapse = "\n")
      
      p <- p + annotate("text",
                        x = 25,
                        y = max(df_filtered$Predicted_LT_ratio, na.rm = TRUE) * 0.8,
                        label = label_text,
                        hjust = 0,
                        size = 3,
                        color = "black")
    }
  }
  
  return(p)
}
p_df1_LT <- plot_turnover_curve(df = df1_Prediction_turnover,df_name = "df1",protein_id = "P11021",model_info_df = result_df1_turnover$Protein_turnover_Rates,raw_points = LT_ratio_filtered) 
print(p_df1_LT)
# ---- End of code block ----



#light, heavy, total AUC  and h/l count(calculate based on predicted AUC)
# ---- Start of code block ----
heavy_count <- heavy_auc %>%
  filter(!is.na(AUC_heavy)) %>%
  group_by(Condition) %>%
  summarise(Protein_Count = n_distinct(PG.ProteinAccessions)) %>%
  mutate(Intensity = "Heavy")

light_count <- light_auc %>%
  filter(!is.na(AUC_light)) %>%
  group_by(Condition) %>%
  summarise(Protein_Count = n_distinct(PG.ProteinAccessions)) %>%
  mutate(Intensity = "Light")

total_count <- total_AUC %>%
  filter(!is.na(AUC_total)) %>%
  group_by(Condition) %>%
  summarise(Protein_Count = n_distinct(PG.ProteinAccessions)) %>%
  mutate(Intensity = "Total")

hl_df <- protein_timepoint_stats %>%
  select(Condition, proteins_1_or_more) %>%
  rename(
    Treatment = Condition,
    Protein_Count = proteins_1_or_more
  ) %>%
  mutate(Intensity = "H/L ratio") %>%
  select(Intensity, Treatment, Protein_Count)

auc_count_all <- bind_rows(heavy_count, light_count, total_count)
auc_count_all <- auc_count_all %>%
  mutate(
    Treatment = recode(Condition,
                       "Tun" = "Tm",
                       "Thg" = "Th",
                       "Con" = "DMSO"),
    Treatment = factor(Treatment, levels = c("DMSO", "Tm", "Th"))
  )
combined_df <- bind_rows(auc_count_all, hl_df)


treat_colors <- c("DMSO" = "#1a80bb", "Tm" = "#ea801c", "Th" = "#b8b8b8")

combined_df$Intensity <- factor(combined_df$Intensity,
                                levels = c("Light", "Heavy","Total", "H/L ratio"))

ggplot(combined_df, aes(x = Treatment, y = Protein_Count, fill = Treatment)) +
  geom_col(width = 0.6, color = "black")+
  geom_text(aes(label = Protein_Count), vjust = -0.3, size = 3.5, color = "black") +
  scale_fill_manual(values = treat_colors) +
  scale_y_continuous(limits = c(0, 10000)) +
  scale_x_discrete(limits = c("DMSO", "Tm", "Th"))+
  facet_wrap(~Intensity, nrow = 1, strip.position = "top") +
  labs(
    x = "Treatment",
    y = "Protein Count"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    theme(panel.spacing = unit(1.5, "cm"))
  )
# ---- End of code block ----



#histogram distribution plot
create_combined_histogram <- function(df1_Prediction_light, df1_Prediction, total_AUC, df1_Prediction_turnover, treatment = "Thg") {

  colors <- list(
    "Thg" = "#b8b8b8",
    "Tun" = "#ea801c"
  )
  fill_color <- colors[[treatment]]
  
  df1_Prediction_light <- df1_Prediction_light %>%
    mutate(log2_AUC = log2(AUC))
  
  df1_Prediction <- df1_Prediction %>%
    mutate(log2_AUC = log2(AUC))
  
  total_AUC <- total_AUC %>%
    mutate(log2_AUC_total = log2(AUC_total))
  
  df1_Prediction_turnover <- df1_Prediction_turnover %>%
    mutate(log2_Turnover_Rate = log2(Turnover_Rate))
  
  # 1. Light (Degradation) data
  df_light <- df1_Prediction_light %>%
    filter(Condition %in% c(treatment, "DMSO"), Time == 0) %>%
    distinct(PG.ProteinAccessions, PG.Genes, Condition, log2_AUC) %>%
    pivot_wider(names_from = Condition, values_from = log2_AUC) %>%
    drop_na(DMSO, all_of(treatment)) %>%
    mutate(log2_ratio = .data[[treatment]] - DMSO,
           Type = "Light (Degradation)",
           n_count = n(),
           median_val = median(log2_ratio, na.rm = TRUE))
  
  # 2. Heavy (Synthesis) data
  df_heavy <- df1_Prediction %>%
    filter(Condition %in% c(treatment, "DMSO"), Time == 0) %>%
    distinct(PG.ProteinAccessions, PG.Genes, Condition, log2_AUC) %>%
    pivot_wider(names_from = Condition, values_from = log2_AUC) %>%
    drop_na(DMSO, all_of(treatment)) %>%
    mutate(log2_ratio = .data[[treatment]] - DMSO,
           Type = "Heavy (Synthesis)",
           n_count = n(),
           median_val = median(log2_ratio, na.rm = TRUE))
  
  # 3. Total data
  df_total <- total_AUC %>%
    filter(Condition %in% c(treatment, "DMSO")) %>%
    distinct(PG.ProteinAccessions, PG.Genes, Condition, log2_AUC_total) %>%
    pivot_wider(names_from = Condition, values_from = log2_AUC_total) %>%
    drop_na(DMSO, all_of(treatment)) %>%
    mutate(log2_ratio = .data[[treatment]] - DMSO,
           Type = "Total",
           n_count = n(),
           median_val = median(log2_ratio, na.rm = TRUE))
  
  # 4. H/L (Turnover rate) data
  df_turnover <- df1_Prediction_turnover %>%
    filter(Treatment %in% c(treatment, "DMSO"), Time == 0) %>%
    distinct(PG.ProteinAccessions, PG.Genes, Treatment, log2_Turnover_Rate) %>%
    pivot_wider(names_from = Treatment, values_from = log2_Turnover_Rate) %>%
    drop_na(DMSO, all_of(treatment)) %>%
    mutate(log2_ratio = .data[[treatment]] - DMSO,
           Type = "H/L (Turnover Rate)",
           n_count = n(),
           median_val = median(log2_ratio, na.rm = TRUE))
  
  combined_data <- bind_rows(
    df_light %>% select(log2_ratio, Type, n_count, median_val),
    df_heavy %>% select(log2_ratio, Type, n_count, median_val),
    df_total %>% select(log2_ratio, Type, n_count, median_val),
    df_turnover %>% select(log2_ratio, Type, n_count, median_val)
  )
  
  combined_data$Type <- factor(combined_data$Type, 
                               levels = c("Light (Degradation)", "Heavy (Synthesis)", 
                                          "H/L (Turnover Rate)","Total"))
  
  annotation_data <- combined_data %>%
    group_by(Type) %>%
    summarise(
      n_count = first(n_count),
      median_val = first(median_val),
      .groups = 'drop'
    )

  p <- ggplot(combined_data, aes(x = log2_ratio)) +
    geom_histogram(binwidth = 0.25, fill = fill_color, color = "black") +
    geom_vline(aes(xintercept = median_val), 
               linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_vline(xintercept = 0, 
               linetype = "dotted", color = "grey40", linewidth = 1) +

    geom_text(data = annotation_data,
              aes(x = 3, y = Inf, 
                  label = paste0("n = ", n_count, "\nmedian = ", round(median_val, 2))),
              vjust = 1, hjust = 1, size = 2) +

    facet_wrap(~Type, nrow = 1, scales = "free_y") +

    scale_x_continuous(limits = c(-3, 3)) +
    scale_y_continuous(limits = c(0, 2500)) +

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
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      strip.text = element_text(size = 10, color = "black"),
      strip.background = element_rect(fill = "white", color = "black")
    )
  
  return(p)
}

p_thg <- create_combined_histogram(df1_Prediction_light, df1_Prediction, total_AUC, df1_Prediction_turnover, treatment = "Thg")
print(p_thg)

p_tun <- create_combined_histogram(df1_Prediction_light, df1_Prediction, total_AUC, df1_Prediction_turnover, treatment = "Tun")
print(p_tun)
