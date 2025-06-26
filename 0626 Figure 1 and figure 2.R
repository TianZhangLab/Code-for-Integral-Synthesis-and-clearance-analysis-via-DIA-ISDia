#0626 Figure 1 and figure 2
library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)
library(VennDiagram)
library(gridExtra)
library(patchwork)  
#DDA data processed
# ---- Start of code block ----
DDA <- read_excel("D:/personal/UVA/Data analysis/202412 SILAC/FIG1. DDA and DIA/DDA_combined_protein_label_quant.xlsx")
DDA <- DDA %>%
  dplyr::select(`Protein ID`,`ratio1_1 Light MaxLFQ Intensity`,`ratio1_2 Light MaxLFQ Intensity`, `ratio1_3 Light MaxLFQ Intensity`,`ratio27_1 Light MaxLFQ Intensity`,
                `ratio27_2 Light MaxLFQ Intensity`, `ratio27_3 Light MaxLFQ Intensity`,`ratio3_1 Light MaxLFQ Intensity`, `ratio3_2 Light MaxLFQ Intensity`,
                `ratio3_3 Light MaxLFQ Intensity`,`ratio9_1 Light MaxLFQ Intensity`,`ratio9_2 Light MaxLFQ Intensity`,`ratio9_3 Light MaxLFQ Intensity`,
                `ratio1_1 Heavy MaxLFQ Intensity`,`ratio1_2 Heavy MaxLFQ Intensity`,`ratio1_3 Heavy MaxLFQ Intensity`,`ratio27_1 Heavy MaxLFQ Intensity`,
                `ratio27_2 Heavy MaxLFQ Intensity`,`ratio27_3 Heavy MaxLFQ Intensity`, `ratio3_1 Heavy MaxLFQ Intensity`, `ratio3_2 Heavy MaxLFQ Intensity`,
                `ratio3_3 Heavy MaxLFQ Intensity`,`ratio9_1 Heavy MaxLFQ Intensity`,`ratio9_2 Heavy MaxLFQ Intensity`,`ratio9_3 Heavy MaxLFQ Intensity`
  )


processed_DDA <- DDA %>%
  pivot_longer(
    cols = -`Protein ID`,  
    names_to = c("condition", "replicate", "type"),
    names_pattern = "ratio(\\d+)_([1-3]) (Light|Heavy) MaxLFQ Intensity",  
    values_to = "Intensity"
  ) %>%
  mutate(
    condition = paste0("ratio", condition),  
    PG_MS2Channel = paste0("PG.MS2Channel", ifelse(type == "Light", 1, 2), replicate) 
  ) %>%
  dplyr::select(`Protein ID`, condition, PG_MS2Channel, Intensity) %>%
  pivot_wider(
    names_from = PG_MS2Channel,  
    values_from = Intensity
  ) %>%
  mutate(Acquisition = "DDA")  %>% 
  dplyr::rename(
    PG.ProteinAccessions = `Protein ID`,  
    R.Condition = condition               
  )



filtered_DDA <- processed_DDA %>%
  rowwise() %>%
  mutate(
    PG.MS2Channel11 = ifelse(is.na(PG.MS2Channel11) | is.na(PG.MS2Channel21), NA, PG.MS2Channel11),
    PG.MS2Channel21 = ifelse(is.na(PG.MS2Channel11) | is.na(PG.MS2Channel21), NA, PG.MS2Channel21),
    PG.MS2Channel12 = ifelse(is.na(PG.MS2Channel12) | is.na(PG.MS2Channel22), NA, PG.MS2Channel12),
    PG.MS2Channel22 = ifelse(is.na(PG.MS2Channel12) | is.na(PG.MS2Channel22), NA, PG.MS2Channel22),
    PG.MS2Channel13 = ifelse(is.na(PG.MS2Channel13) | is.na(PG.MS2Channel23), NA, PG.MS2Channel13),
    PG.MS2Channel23 = ifelse(is.na(PG.MS2Channel13) | is.na(PG.MS2Channel23), NA, PG.MS2Channel23)
  ) %>%
  filter(
    !(is.na(PG.MS2Channel11) & is.na(PG.MS2Channel12) & is.na(PG.MS2Channel13) &
        is.na(PG.MS2Channel21) & is.na(PG.MS2Channel22) & is.na(PG.MS2Channel23))
  )%>%
  ungroup()
# ---- End of code block ----


#DIA data process
# ---- Start of code block ----
clean_data <- function(df) {
  df %>%
    mutate(
      PG.MS2Channel1 = na_if(PG.MS2Channel1, "NaN"),
      PG.MS2Channel2 = na_if(PG.MS2Channel2, "NaN")
    ) %>%
    mutate(
      PG.MS2Channel1 = ifelse(is.na(PG.MS2Channel1) | is.na(PG.MS2Channel2), NA, PG.MS2Channel1),
      PG.MS2Channel2 = ifelse(is.na(PG.MS2Channel1) | is.na(PG.MS2Channel2), NA, PG.MS2Channel2)
    ) %>%
    filter(!is.na(PG.MS2Channel1) & !is.na(PG.MS2Channel2))
}

DIA_noFAIMS_win2_IT3_5 <- read_excel("D:/personal/UVA/Data analysis/202412 SILAC/FIG1. DDA and DIA/DIA_noFAIMS_win2_IT3.5.xlsx")
cleaned_DIA <- clean_data(DIA_noFAIMS_win2_IT3_5)


#reshape df
reshaped_DIA <- cleaned_DIA %>%
  mutate(R.FileName = case_when(
    R.FileName %in% c("156_27_1_noFAIMS_DIA", "153_3_1_noFAIMS_DIA","026_9_1","038_1_1") ~ "1",
    R.FileName %in% c("157_27_2_noFAIMS_DIA", "154_3_2_noFAIMS_DIA","027_9_2","039_1_2_20241220083228") ~ "2",
    R.FileName %in% c("158_27_3_noFAIMS_DIA", "155_3_3_noFAIMS_DIA","028_9_3","040_1_3") ~ "3",
    TRUE ~ as.character(R.FileName)  
  ))%>%
  pivot_wider(
    names_from = R.FileName,
    values_from = c(PG.MS2Channel1, PG.MS2Channel2),
    names_glue = "{.value}{R.FileName}"  
  )  %>%
  mutate(Acquisition = "DIA")%>%
  filter(
    !(is.na(PG.MS2Channel11) & is.na(PG.MS2Channel12) & is.na(PG.MS2Channel13) &
        is.na(PG.MS2Channel21) & is.na(PG.MS2Channel22) & is.na(PG.MS2Channel23))
  )
# ---- End of code block ----


# combine DDA and DIA
columns_to_convert <- c("PG.MS2Channel11", "PG.MS2Channel12", "PG.MS2Channel13",
                        "PG.MS2Channel21", "PG.MS2Channel22", "PG.MS2Channel23")
filtered_DDA[columns_to_convert] <- lapply(filtered_DDA[columns_to_convert], as.numeric)
reshaped_DIA[columns_to_convert] <- lapply(reshaped_DIA[columns_to_convert], as.numeric)

combined <- bind_rows(filtered_DDA, reshaped_DIA)

# calculate CV  and ratio show in 3 replicates
# ---- Start of code block ----
calculate_CV_and_LH <- function(df) {
  df <- df %>%
    # Calculate light intensity CV
    rowwise() %>%
    mutate(
      Light_mean = mean(c(`PG.MS2Channel11`, `PG.MS2Channel12`, `PG.MS2Channel13`), na.rm = TRUE),
      CV_Light = {
        values <- c_across(c(`PG.MS2Channel11`, `PG.MS2Channel12`, `PG.MS2Channel13`))
        values <- values[!is.na(values)]
        if (length(values) >= 3) {
          (sd(values) / mean(values)) * 100
        } else {
          NA
        }
      }
    ) %>%
    ungroup() %>%
    
    # Calculate heavy intensity CV
    rowwise() %>%
    mutate(
      Heavy_mean = mean(c(`PG.MS2Channel21`, `PG.MS2Channel22`, `PG.MS2Channel23`), na.rm = TRUE),
      CV_Heavy = {
        values <- c_across(c(`PG.MS2Channel21`, `PG.MS2Channel22`, `PG.MS2Channel23`))
        values <- values[!is.na(values)]
        if (length(values) >= 3) {
          (sd(values) / mean(values)) * 100
        } else {
          NA
        }
      }
    ) %>%
    ungroup() %>%
    
    # Calculate total intensity and CV
    mutate(
      Total_intensity1 = PG.MS2Channel11 + PG.MS2Channel21,
      Total_intensity2 = PG.MS2Channel12 + PG.MS2Channel22,
      Total_intensity3 = PG.MS2Channel13 + PG.MS2Channel23
    ) %>%
    rowwise() %>%
    mutate(
      Total_intensity_CV = {
        values <- c(Total_intensity1, Total_intensity2, Total_intensity3)
        values <- values[!is.na(values)]
        if (length(values) >= 3) {
          (sd(values) / mean(values)) * 100
        } else {
          NA_real_
        }
      }
    ) %>%
    ungroup() %>%
    
    # Calculate H/L ratio mean and CV
    mutate(
      `H/L ratio 1` = `PG.MS2Channel21` / `PG.MS2Channel11`,
      `H/L ratio 2` = `PG.MS2Channel22` / `PG.MS2Channel12`,
      `H/L ratio 3` = `PG.MS2Channel23` / `PG.MS2Channel13`
    ) %>%
    rowwise() %>%
    mutate(
      HtoL_ratio_mean = mean(c(`H/L ratio 1`, `H/L ratio 2`, `H/L ratio 3`), na.rm = TRUE),
      HtoL_ratio_CV = {
        values <- c(`H/L ratio 1`, `H/L ratio 2`, `H/L ratio 3`)
        values <- values[!is.na(values)]
        if (length(values) >= 3) {
          (sd(values) / mean(values)) * 100
        } else {
          NA_real_
        }
      }
    )%>%
    ungroup()
  
  return(df)
}
calculated_combined <- calculate_CV_and_LH(combined)
# ---- End of code block ----

#Figure 1
# count PG numbers between DDA and DIA at different H/L ratio
# ---- Start of code block ----
columns_to_count <- c("PG.MS2Channel11", "PG.MS2Channel12", "PG.MS2Channel13", 
                      "PG.MS2Channel21", "PG.MS2Channel22", "PG.MS2Channel23")

result <- calculated_combined %>%
  filter(R.Condition %in% c("ratio1", "ratio3", "ratio9", "ratio27")) %>%
  rowwise() %>%
  mutate(num_numeric = sum(!is.na(c_across(all_of(columns_to_count))))) %>%
  group_by(Acquisition, R.Condition) %>%
  summarise(
    show_in_1_replicates = sum(num_numeric >= 2),
    show_in_2_replicates = sum(num_numeric >= 4),
    show_in_3_replicates = sum(num_numeric >= 6),
    # Counting PG.ProteinAccessions under different CV conditions
    LightCV_under_20_percent = n_distinct(PG.ProteinAccessions[
      num_numeric == 6 & CV_Light < 20 & !is.na(HtoL_ratio_CV)
    ]),
    HeavyCV_under_20_percent = n_distinct(PG.ProteinAccessions[
      num_numeric == 6 & CV_Heavy < 20 & !is.na(HtoL_ratio_CV)
    ]),
    Total_intensity_CV_under_20_percent = n_distinct(PG.ProteinAccessions[
      num_numeric == 6 & Total_intensity_CV < 20 & !is.na(HtoL_ratio_CV)
    ]),
    Intensity_CV_under_20_percent = n_distinct(PG.ProteinAccessions[
      num_numeric == 6 & CV_Light < 20 &CV_Heavy < 20 & Total_intensity_CV < 20 & HtoL_ratio_CV<20
    ]),
    .groups = "drop"
  )


# Define the order of categories
stack_columns <- c("Intensity_CV_under_20_percent", "show_in_3_replicates", 
                   "show_in_2_replicates", "show_in_1_replicates")


# Compute the incremental stacking differences
long_data <- result %>%
  mutate(
    bottom = Intensity_CV_under_20_percent,  
    second = show_in_3_replicates - bottom,  
    third = show_in_2_replicates - show_in_3_replicates,  
    top = show_in_1_replicates - show_in_2_replicates  
  ) %>%
  pivot_longer(
    cols = c(bottom, second, third, top),
    names_to = "Category",
    values_to = "Incremental_Count"
  ) %>%
  mutate(
    Category = factor(Category, 
                      levels = c("bottom", "second", "third", "top"),
                      labels = c("Intensity_CV_under_20_percent", 
                                 "show_in_3_replicates", 
                                 "show_in_2_replicates", 
                                 "show_in_1_replicates"))
  ) %>%
  group_by(R.Condition, Acquisition) %>%
  mutate(
    label = cumsum(Incremental_Count)
  ) %>%
  ungroup() %>%
  mutate(
    Category = factor(Category, levels = rev(levels(Category))),
    R.Condition = factor(R.Condition, levels = c("ratio1", "ratio3", "ratio9", "ratio27"))
  )


library(colorspace)

ratio_hues <- c(
  "ratio1" = "#00cd6c",  # green
  "ratio3" = "#009ade",  # blue
  "ratio9" = "#af58ba",  # purple
  "ratio27" = "#ffc61e"   # yellow
)


category_order <- c(
  "Intensity_CV_under_20_percent",
  "show_in_3_replicates",
  "show_in_2_replicates",        
  "show_in_1_replicates"
)

create_color_palette <- function(ratio_hues, category_order) {
  palette <- c()
  
  brightness_adjust <- c(0, 0.35, 0.65, 0.95)  
  
  for (ratio_name in names(ratio_hues)) {
    base_color <- ratio_hues[ratio_name]
    
    for (i in 1:length(category_order)) {
      category <- category_order[i]
      key <- paste(ratio_name, category, sep = "_")
      
      if (i == 1) {
        # show_in_1_replicates 
        palette[key] <- base_color
      } else {

        palette[key] <- lighten(base_color, amount = brightness_adjust[i])
      }
    }
  }
  
  return(palette)
}


ggplot(long_data, aes(x = R.Condition, y = Incremental_Count, fill = interaction(R.Condition, Category, sep = "_"))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) + 
  facet_wrap(~ Acquisition) + 
  scale_x_discrete(labels = c("ratio1" = "1/1", "ratio3" = "3/1", "ratio9" = "9/1", "ratio27" = "27/1")) +
  geom_text(data = subset(long_data, Acquisition=="DIA" & Category == "show_in_1_replicates"),
            aes(label = label),
            position = position_nudge(y = 10300),  
            size = 3, color = "black")+
  geom_text(aes(label = ifelse(Acquisition == "DIA" & Category == "show_in_1_replicates", NA, label)),
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black") + 
  scale_fill_manual(guide = "none",
                    values = create_color_palette(ratio_hues, category_order),
                    labels = function(breaks) {
                      parts <- strsplit(breaks, "_")
                      ratio_names <- sapply(parts, function(x) x[1])
                      ratio_labels <- c("ratio1" = "1/1", "ratio3" = "3/1", "ratio9" = "9/1", "ratio27" = "27/1")
                      
                      category_names <- sapply(parts, function(x) {
                        paste(x[-(1:1)], collapse = "_")  
                      })
                      category_labels <- c("Intensity_CV_under_20_percent" = "3 (CV<20%)",
                                           "show_in_3_replicates" = "3",
                                           "show_in_2_replicates" = "2",
                                           "show_in_1_replicates" = "1")
                      
                      paste(ratio_labels[ratio_names], "-", category_labels[category_names])
                    }
  ) +
  scale_y_continuous(breaks = c(0,5000, 10000), limits = c(0, 11000)) +  
  labs(
    title = NULL,
    x = "H/L Ratio",
    y = "Protein count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black",hjust = 0.5, vjust = 1,size = 12),  
    axis.text.y = element_text(color = "black",size = 12),  
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14), 
    strip.text = element_text(size = 12,face = "bold"),  
    legend.text = element_text(size = 10) , 
    legend.position = c(0.1, 0.9), 
    legend.justification = c(0, 1) , 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", size = 0.5),  
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.3, "cm")
  ) +
  geom_hline(yintercept = 0, color = "white", size = 0.5) + 
  geom_vline(xintercept = 0, color = "black", size = 1) 
# ---- End of code block ----

#Light, heavy, total and ratio CV comparison
# ---- Start of code block ----
long_combined <- calculated_combined %>%
  dplyr::select(R.Condition, Acquisition,CV_Light, CV_Heavy, Total_intensity_CV, HtoL_ratio_CV) %>%
  pivot_longer(cols = c(CV_Light, CV_Heavy, Total_intensity_CV, HtoL_ratio_CV),
               names_to = "Group", values_to = "CV")
long_combined <- long_combined %>%
  mutate(R.Condition = factor(R.Condition, levels = c("ratio1", "ratio3", "ratio9", "ratio27")))

long_combined <- long_combined %>%
  mutate(Group = factor(Group, levels = c("CV_Light", "CV_Heavy","Total_intensity_CV", "HtoL_ratio_CV")))

group_labels <- c(
  "Total_intensity_CV" = "Total Intensity",
  "HtoL_ratio_CV" = "H/L Ratio",
  "CV_Light" = "Light",
  "CV_Heavy" = "Heavy"
)

# boxplot of different CV at different ratio
long_combined$Group <- as.factor(long_combined$Group)

group_list <- unique(long_combined$Group)

plot_list_1 <- list()
for (g in group_list) {
  data_subset <- long_combined %>% filter(Group == g)
  
  title_name <- case_when(
    g == "CV_Light" ~ "Light intensity",
    g == "CV_Heavy" ~ "Heavy intensity",
    g == "HtoL_ratio_CV" ~ "H/L ratio",
    g == "Total_intensity_CV" ~ "Total intensity",
    TRUE ~ as.character(g) 
  )
  
  signif_data <- data_subset %>%
    group_by(R.Condition) %>%
    summarise(
      p_value = tryCatch(
        t.test(CV ~ Acquisition)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        is.na(p_value)     ~ "NA",
        p_value <= 0.001   ~ "***",
        p_value <= 0.01    ~ "**",
        p_value <= 0.05    ~ "*",
        TRUE               ~ "ns"
      )
    )
  
  p <- ggplot(data_subset, aes(x = R.Condition, y = CV, fill = as.factor(Acquisition))) +
    geom_boxplot(width = 0.75, color = "black", outlier.shape = NA,
                 position = position_dodge(width = 0.75)) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 18, size = 3,
                 position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = 20, color = "red", linetype = "dotted", size = 1) +
    geom_text(
      data = signif_data,
      aes(x = R.Condition, y = 80, label = label),
      inherit.aes = FALSE,
      size = 5,
      color = "black"
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_x_discrete(
      labels = c("ratio1" = "1/1", "ratio3" = "3/1", "ratio9" = "9/1", "ratio27" = "27/1"),
      expand = expansion(mult = c(0.3, 0.3))
    ) +
    scale_fill_manual(
      values = c("DDA" = "#Ff1f5b", "DIA" = "#009ade"),
      labels = c("DDA" = "DDA", "DIA" = "DIA")
    ) +
    labs(
      title = title_name,
      x = "H/L ratio",
      y = "Protein group (CV%)",
      fill = "Acquisition"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(color = "black", hjust = 0.5, vjust = 1, size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 14,  color = "black"),
      axis.title.y = element_text(size = 14,  color = "black"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 12, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      panel.border = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.3, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
    )
  
  plot_list_1[[g]] <- p
}

wrap_plots(plot_list_1, nrow = 1)
# ---- End of code block ----

#Venn diagram and histogram plot of intensity to understand the difference between DIA and DDA
# ---- Start of code block ----
columns_to_count <- c(
  "PG.MS2Channel11", "PG.MS2Channel12", "PG.MS2Channel13",
  "PG.MS2Channel21", "PG.MS2Channel22", "PG.MS2Channel23"
)

filtered_combined <- calculated_combined %>%
  rowwise() %>%
  mutate(num_numeric = sum(!is.na(c_across(any_of(columns_to_count)))))  %>%
  ungroup() %>%
  filter(
    num_numeric == 6,
    CV_Light < 20,
    CV_Heavy < 20,
    Total_intensity_CV < 20,
    HtoL_ratio_CV<20
  ) %>%
  distinct(PG.ProteinAccessions, R.Condition, Acquisition,PG.MS2Channel11,PG.MS2Channel12,PG.MS2Channel13,PG.MS2Channel21,PG.MS2Channel22,PG.MS2Channel23,CV_Light,CV_Heavy,Total_intensity_CV)%>%
  mutate(
    Light_mean = rowMeans(
      across(c("PG.MS2Channel11", "PG.MS2Channel12", "PG.MS2Channel13")),
      na.rm = TRUE
    ),
    Heavy_mean = rowMeans(
      across(c("PG.MS2Channel21", "PG.MS2Channel22", "PG.MS2Channel23")),
      na.rm = TRUE
    ),
    Log10_Light_mean = log10(Light_mean) ,
    Log10_Heavy_mean =log10(Heavy_mean)
  )

set1f_H1L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio1"& filtered_combined$Acquisition == "DIA"])
set2f_H1L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio1"& filtered_combined$Acquisition == "DDA"])
set1f_H3L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio3"& filtered_combined$Acquisition == "DIA"])
set2f_H3L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio3"& filtered_combined$Acquisition == "DDA"])
set1f_H9L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio9"& filtered_combined$Acquisition == "DIA"])
set2f_H9L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio9"& filtered_combined$Acquisition == "DDA"])
set1f_H27L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio27"& filtered_combined$Acquisition == "DIA"])
set2f_H27L1 <- unique(filtered_combined$PG.ProteinAccessions[filtered_combined$R.Condition == "ratio27"& filtered_combined$Acquisition == "DDA"])

#Venn plot
# H1L1 Venn plot
grid.newpage()
grid.text("H/L = 1/1", x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))
size1 <- length(set1f_H1L1)  
size2 <- length(set2f_H1L1)   
overlap <- length(intersect(set1f_H1L1, set2f_H1L1))
circle1_size <- sqrt(size1) 
circle2_size <- sqrt(size2)

# Normalize to make the larger circle bigger
scale_factor <- max(circle1_size, circle2_size) / 0.4  
circle1_size <- circle1_size / scale_factor
circle2_size <- circle2_size*1.8 / scale_factor

grid.circle(x = 0.4, y = 0.5, r = circle1_size, gp = gpar(fill = "#009ade", alpha = 0.7, col = "black", lwd = 3))
grid.circle(x = 0.6, y = 0.5, r = circle2_size, gp = gpar(fill = "#Ff1f5b", alpha = 0.7, col = "black", lwd = 3))

grid.text(size1 - overlap, x = 0.15, y = 0.5, gp = gpar(fontsize = 18))  # DIA 专属蛋白
grid.text(size2 - overlap, x = 0.87, y = 0.5, gp = gpar(fontsize = 18))  # DDA 专属蛋白
grid.text(overlap, x = 0.5, y = 0.5, gp = gpar(fontsize = 18))  # 重叠蛋白

grid.text("DIA", x = 0.1, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("DDA", x = 0.75, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))


#H3L1 Venn plot
grid.newpage()
grid.text("H/L = 3/1", x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))
size1 <- length(set1f_H3L1)   
size2 <- length(set2f_H3L1)  
overlap <- length(intersect(set1f_H3L1, set2f_H3L1))

circle1_size <- sqrt(size1)  
circle2_size <- sqrt(size2)


scale_factor <- max(circle1_size, circle2_size) / 0.4  
circle1_size <- circle1_size / scale_factor
circle2_size <- circle2_size*1.8 / scale_factor

grid.circle(x = 0.4, y = 0.5, r = circle1_size, gp = gpar(fill = "#009ade", alpha = 0.7, col = "black", lwd = 3))
grid.circle(x = 0.6, y = 0.5, r = circle2_size, gp = gpar(fill = "#Ff1f5b", alpha = 0.7, col = "black", lwd = 3))

grid.text(size1 - overlap, x = 0.15, y = 0.5, gp = gpar(fontsize = 18)) 
grid.text(size2 - overlap, x = 0.87, y = 0.5, gp = gpar(fontsize = 18))  
grid.text(overlap, x = 0.5, y = 0.5, gp = gpar(fontsize = 18))  

grid.text("DIA", x = 0.1, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("DDA", x = 0.75, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))

#H9L1 Venn plot
grid.newpage()
grid.text("H/L = 9/1", x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))
size1 <- length(set1f_H9L1)  
size2 <- length(set2f_H9L1)   
overlap <- length(intersect(set1f_H9L1, set2f_H9L1))

circle1_size <- sqrt(size1)  
circle2_size <- sqrt(size2)

scale_factor <- max(circle1_size, circle2_size) / 0.4  
circle1_size <- circle1_size / scale_factor
circle2_size <- circle2_size*1.9 / scale_factor

grid.circle(x = 0.4, y = 0.5, r = circle1_size, gp = gpar(fill = "#009ade", alpha = 0.7, col = "black", lwd = 3))
grid.circle(x = 0.6, y = 0.5, r = circle2_size, gp = gpar(fill = "#Ff1f5b", alpha = 0.7, col = "black", lwd = 3))

grid.text(size1 - overlap, x = 0.15, y = 0.5, gp = gpar(fontsize = 18))  
grid.text(size2 - overlap, x = 0.8, y = 0.5, gp = gpar(fontsize = 18))  
grid.text(overlap, x = 0.6, y = 0.5, gp = gpar(fontsize = 18))  


grid.text("DIA", x = 0.1, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("DDA", x = 0.75, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))

#H27L1 Venn plot
grid.newpage()
grid.text("H/L = 27/1", x = 0.5, y = 0.95, gp = gpar(fontsize = 18, fontface = "bold"))
size1 <- length(set1f_H27L1)  
size2 <- length(set2f_H27L1)   
overlap <- length(intersect(set1f_H27L1, set2f_H27L1))

circle1_size <- sqrt(size1)  
circle2_size <- sqrt(size2)

scale_factor <- max(circle1_size, circle2_size) / 0.4  
circle1_size <- circle1_size / scale_factor
circle2_size <- circle2_size*1.9 / scale_factor

grid.circle(x = 0.4, y = 0.5, r = circle1_size, gp = gpar(fill = "#009ade", alpha = 0.7, col = "black", lwd = 3))
grid.circle(x = 0.6, y = 0.5, r = circle2_size, gp = gpar(fill = "#Ff1f5b", alpha = 0.7, col = "black", lwd = 3))

grid.text(size1 - overlap, x = 0.15, y = 0.5, gp = gpar(fontsize = 18)) 
grid.text(size2 - overlap, x = 0.8, y = 0.5, gp = gpar(fontsize = 18))  
grid.text(overlap, x = 0.6, y = 0.5, gp = gpar(fontsize = 18))  

grid.text("DIA", x = 0.1, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("DDA", x = 0.75, y = 0.9, gp = gpar(fontsize = 14, fontface = "bold"))

# ---- End of code block ----


#histogram overlap and non overlap
common_proteins_ratio1 <- intersect(set1f_H1L1 , set2f_H1L1 )
common_proteins_ratio3 <- intersect(set1f_H3L1 , set2f_H3L1 )
common_proteins_ratio9 <- intersect(set1f_H9L1 , set2f_H9L1 )
common_proteins_ratio27 <- intersect(set1f_H27L1 , set2f_H27L1 )

#Abundance distribution
#DIA light
# ---- Start of code block ----
binned_data_2 <- filtered_combined %>%
  filter(Acquisition == "DIA", !is.na(Log10_Light_mean)) %>%
  mutate(
    bin = cut(
      Log10_Light_mean,
      breaks = c(-Inf, seq(1, 6.8, by = 0.2), Inf),
      include.lowest = TRUE
    )
  ) %>%
  group_by(bin) %>%
  summarise(
    Count = n(),
    Proteins = list(unique(`PG.ProteinAccessions`)),
    .groups = "drop"
  )

label_overlap_2 <- filtered_combined %>%
  filter(Acquisition == "DIA", !is.na(Log10_Light_mean)) %>%
  mutate(
    bin = cut(
      Log10_Light_mean,
      breaks = c(-Inf, seq(1, 6.8, by = 0.2), Inf),
      include.lowest = TRUE
    ),
    DDA_overlap = case_when(
      R.Condition == "ratio1" & PG.ProteinAccessions %in% common_proteins_ratio1 ~ "overlap",
      R.Condition == "ratio3" & PG.ProteinAccessions %in% common_proteins_ratio3 ~ "overlap",
      R.Condition == "ratio9" & PG.ProteinAccessions %in% common_proteins_ratio9 ~ "overlap",
      R.Condition == "ratio27" & PG.ProteinAccessions %in% common_proteins_ratio27 ~ "overlap",
      TRUE ~ "non overlap"
    )
  )

bin_summary_2 <- label_overlap_2 %>%
  group_by(bin, R.Condition) %>%
  mutate(Total_Count = n()) %>%
  group_by(bin, R.Condition, DDA_overlap) %>%
  summarise(
    Count       = n(),
    Total_Count = Total_Count[1], 
    Percentage  = (Count / Total_Count[1]) * 100, 
    .groups     = "drop"
  ) %>%
  mutate(
    Label = case_when(
      bin == "(1,1.2]"  ~ "1",
      bin == "(2,2.2]"  ~ "2",
      bin == "(3,3.2]"  ~ "3",
      bin == "(4,4.2]"  ~ "4",
      bin == "(5,5.2]"  ~ "5",
      bin == "(6,6.2]"  ~ "6",
      TRUE             ~ NA_character_
    )
  )

bin_summary_2 <- bin_summary_2 %>%
  mutate(bin = as.character(bin))

lookup <- bin_summary_2 %>%
  filter(!is.na(Label)) %>%     
  distinct(bin, Label)          

bin_summary_2 <- bin_summary_2 %>%
  mutate(R.Condition = factor(R.Condition, levels = c("ratio1", "ratio3", "ratio9", "ratio27")))

#median
median_real <- label_overlap_2 %>%
  group_by(R.Condition, DDA_overlap) %>%
  summarise(
    median_log10 = median(Log10_Light_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y = ifelse(DDA_overlap == "overlap", 1200, 1150),
    x = Inf,
    hjust = 1.1,
    vjust = 1.5,
    label = paste0("median (", DDA_overlap, ") = ", round(median_log10, 2))
  )%>%
  mutate(R.Condition = factor(R.Condition, levels = c("ratio1", "ratio3", "ratio9", "ratio27")))%>%
  mutate(
    label = case_when(
      DDA_overlap == "non overlap" ~ paste0("Median (DIA only) = ", round(median_log10, 2)),
      DDA_overlap == "overlap"     ~ paste0("Median (Shared) = ", round(median_log10, 2)),
      TRUE                         ~ label
    )
  )


ggplot(bin_summary_2, aes(x = bin, y = Count, fill = DDA_overlap)) +
  geom_col(position = "stack", width = 0.9) +
  scale_y_continuous(
    breaks = c(300, 600, 900),
    limits = c(-50, 1200),   
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    breaks = lookup$bin,    
    labels = lookup$Label   
  ) +
  geom_text(
    data = median_real,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2,
    hjust = median_real$hjust,
    vjust = median_real$vjust
  ) +
  labs(x="Log10(light intensity)", y="Protein count", fill="DIA Overlap") +
  facet_wrap(~ R.Condition, ncol = 4) +  
  scale_fill_manual(
    values = c("non overlap" = "#009ade", "overlap" = "#C94D86"),
    labels = c("non overlap" = "DIA only", "overlap" = "Shared")
  )+
  theme_minimal() +
  theme(
    panel.border       = element_rect(color = "black", fill = NA, size = 0.75),
    panel.background   = element_blank(),
    panel.grid         = element_blank(),
    axis.line.y        = element_line(color = "black"),
    axis.ticks.y       = element_line(color = "black"),
    axis.ticks.x       = element_line(color = "black"),
    axis.title.y       = element_text(color = "black"),
    plot.margin        = margin(t = 20, r = 20, b = 40, l = 20)
  ) +
  coord_cartesian(clip = "off") 
# ---- End of code block ----

#Figure 2
#scatter and density plots
#Before filtering with CV<20%
plot_scatter_with_density <- function(df,
                                      acquisition = "DDA",
                                      intensity_type = c("Light", "Heavy"),
                                      intensity_mean_col = NULL,
                                      color_palette = NULL) {

  intensity_type <- match.arg(intensity_type)
  
  
  if (is.null(intensity_mean_col)) {
    intensity_mean_col <- if (intensity_type == "Light") "Light_mean" else "Heavy_mean"
  }
  
  
  if (is.null(color_palette)) {
    color_palette <- c("1" = "#00CD6C", "3" = "#009ADE", "9" = "#AF58BA", "27" = "#FFC61E")
  }
  
  
  processed_df <- df %>%
    mutate(
      log2_HL = log2(HtoL_ratio_mean),
      log2_Intensity = log2(.data[[intensity_mean_col]])
    ) %>%
    select(PG.ProteinAccessions, R.Condition, Acquisition,
           log2_HL, log2_Intensity,
           CV_Light, CV_Heavy, Total_intensity_CV, HtoL_ratio_CV)
  
  
  filtered_data <- processed_df %>%
    filter(Acquisition == acquisition) %>%
    filter(log2_HL >= -4 & log2_HL <= 6.5) %>%
    mutate(R.Condition = factor(
      as.character(R.Condition),
      levels = c("ratio1", "ratio3", "ratio9", "ratio27"),
      labels = c("1", "3", "9", "27")
    ))
  
  
  log2_lines <- data.frame(y = log2(c(1, 3, 9, 27)))
  
  
  main_plot <- ggplot(filtered_data, aes(x = log2_Intensity, y = log2_HL, color = R.Condition)) +
    geom_point(alpha = 0.5, size = 0.2) +
    coord_cartesian(ylim = c(-4, 6.5)) +
    scale_color_manual(values = color_palette) +
    geom_hline(data = log2_lines, aes(yintercept = y),
               linetype = "dashed", color = "black", size = 0.7) +
    labs(
      x = paste0("log2(", intensity_type, " intensity)"),
      y = "log2(H/L ratio)",
      color = "Condition"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
  
  
  sd_data <- filtered_data %>%
    group_by(R.Condition) %>%
    summarise(
      sd_value = sd(log2_HL, na.rm = TRUE),
      mean_log2_HL = mean(log2_HL, na.rm = TRUE),
      .groups = "drop"
    )
  
  max_x <- max(density(filtered_data$log2_HL, na.rm = TRUE)$y)
  x_pos <- max_x * 5.5
  
  
  density_plot <- ggplot(filtered_data, aes(y = log2_HL, color = R.Condition)) +
    geom_density(size = 1, fill = NA) +
    scale_color_manual(values = color_palette) +
    geom_text(
    data = sd_data,
    aes(x = x_pos, y = mean_log2_HL, label = round(sd_value, 2)),
    inherit.aes = FALSE,
    color = "black",
    size = 3,
    hjust = 0
    ) +
    annotate("text", x = x_pos, y = max(sd_data$mean_log2_HL) + 1,
    label = "s.d.", hjust = 0, size = 4, color = "black") +
    geom_hline(data = log2_lines, aes(yintercept = y),
               linetype = "dashed", color = "gray", size = 0.7) +
    labs(x = "Density", y = NULL, color = "Condition") +
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
    coord_cartesian(xlim = c(0, x_pos * 1.3))
  
  
  final_plot <- main_plot + density_plot +
    patchwork::plot_layout(widths = c(1, 1)) +
    patchwork::plot_annotation(
      title = acquisition,
      theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
    )
  
  return(final_plot)
}


# DDA 
plot_scatter_with_density(calculated_combined, acquisition = "DDA", intensity_type = "Light")

# DIA 
plot_scatter_with_density(calculated_combined, acquisition = "DIA", intensity_type = "Light")

#After filtering with CV<20%
processed_combined_cv20 <-calculated_combined%>%
  filter(CV_Light<20,CV_Heavy<20,Total_intensity_CV<20,HtoL_ratio_CV<20)
# DDA 
plot_scatter_with_density(processed_combined_cv20, acquisition = "DDA", intensity_type = "Light")

# DIA 
plot_scatter_with_density(processed_combined_cv20, acquisition = "DIA", intensity_type = "Light")

#Intensity comparison between different ratio
# ---- Start of code block ----
processed_data_CV_filtered <- calculated_combined %>%
  pivot_longer(
    cols = starts_with("PG.MS2Channel"),  
    names_to = c("Type", "Replicate"),
    names_pattern = "PG\\.MS2Channel(\\d)(\\d+)",  
    values_to = "Intensity"
  ) %>%
  mutate(
    Intensity_Type = ifelse(Type == "1", "Light", "Heavy")
  )  %>%
  filter(CV_Light < 20, CV_Heavy < 20, Total_intensity_CV < 20,HtoL_ratio_CV<20)  

# Compute mean intensity per condition and type
mean_intensity_data_CVfiltered <- processed_data_CV_filtered %>%
  filter(Acquisition == "DIA") %>%
  group_by(R.Condition, Intensity_Type, PG.ProteinAccessions) %>%
  summarise(mean_intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Extract reference values for ratio1
reference_values_CVfiltered <- mean_intensity_data_CVfiltered %>%
  filter(R.Condition == "ratio1") %>%
  select(Intensity_Type, PG.ProteinAccessions, reference_mean = mean_intensity)

# Join reference values back to the main dataset
processed_data_DIA_CVfiltered <- mean_intensity_data_CVfiltered %>%
  left_join(reference_values_CVfiltered, by = c("Intensity_Type", "PG.ProteinAccessions")) %>%
  mutate(Normalized_intensity = mean_intensity / reference_mean)%>%
  mutate(R.Condition = factor(R.Condition, levels = c("ratio1", "ratio3", "ratio9", "ratio27"))) 


ggplot(processed_data_DIA_CVfiltered, aes(x = R.Condition, y = log2(Normalized_intensity), fill = Intensity_Type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +  
  scale_x_discrete(labels = c("ratio1" = "1/1", "ratio3" = "3/1", "ratio9" = "9/1", "ratio27" = "27/1")) +
  scale_fill_manual(values = c("Light" = "#404041", "Heavy" = "#EE1C25")) +  
  stat_summary(
    fun = mean, geom = "point", shape = 18, size = 3, color = "black", 
    position = position_dodge(width = 0.9)  
  ) +
  stat_summary(
    fun = mean, geom = "text", aes(label = round(..y.., 2)), 
    position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5  
  ) +
  labs(
    title = "Comparison of Intensity by Condition in DIA filtered with CV<20%",
    x = "H/L ratio",
    y = "log2(Relative abundance)",
    fill = "Intensity Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 1,size = 10),  
    axis.text.y = element_text(size = 10),  
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10), 
    legend.text = element_text(size = 10) , 
    legend.title = element_blank(),  
    legend.position = c(0.1, 0.1),  
    legend.justification = c(0, 0) ,  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", size = 0.5),  
    panel.border = element_blank(),  
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.3, "cm")  
  ) +
  geom_vline(xintercept = 0, color = "black", size = 1) + 
  coord_cartesian(ylim = c(-7, 3))  

#DDA plot
# Compute mean intensity per condition and type
mean_intensity_data_DDA_CVfiltered <- processed_data_CV_filtered %>%
  filter(Acquisition == "DDA") %>%
  group_by(R.Condition, Intensity_Type, PG.ProteinAccessions) %>%
  summarise(mean_intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Extract reference values for ratio1
reference_values_CVfiltered <- mean_intensity_data_DDA_CVfiltered %>%
  filter(R.Condition == "ratio1") %>%
  select(Intensity_Type, PG.ProteinAccessions, reference_mean = mean_intensity)

# Join reference values back to the main dataset
processed_data_DDA_CVfiltered <- mean_intensity_data_DDA_CVfiltered %>%
  left_join(reference_values_CVfiltered, by = c("Intensity_Type", "PG.ProteinAccessions")) %>%
  mutate(Normalized_intensity = mean_intensity / reference_mean)%>%
  mutate(R.Condition = factor(R.Condition, levels = c("ratio1", "ratio3", "ratio9", "ratio27"))) 



ggplot(processed_data_DDA_CVfiltered, aes(x = R.Condition, y = log2(Normalized_intensity), fill = Intensity_Type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9)) +  
  scale_x_discrete(labels = c("ratio1" = "1/1", "ratio3" = "3/1", "ratio9" = "9/1", "ratio27" = "27/1")) +
  scale_fill_manual(values = c("Light" = "#404041", "Heavy" = "#EE1C25")) +  
  stat_summary(
    fun = mean, geom = "point", shape = 18, size = 3, color = "black", 
    position = position_dodge(width = 0.9)  
  ) +
  stat_summary(
    fun = mean, geom = "text", aes(label = round(..y.., 2)), 
    position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5  
  ) +
  labs(
    title = "Comparison of Intensity by Condition in DDA filtered with CV<20%",
    x = "H/L ratio",
    y = "log2(Relative abundance)",
    fill = "Intensity Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 1,size = 10),  
    axis.text.y = element_text(size = 10),  
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10), 
    legend.text = element_text(size = 10) , 
    legend.title = element_blank(),  
    legend.position = c(0.1, 0.1),  
    legend.justification = c(0, 0) ,  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", size = 0.5),  
    panel.border = element_blank(),  
    axis.ticks = element_line(color = "black", size = 0.5), 
    axis.ticks.length = unit(0.3, "cm")  
  ) +
  geom_vline(xintercept = 0, color = "black", size = 1) + 
  coord_cartesian(ylim = c(-7, 3))  
# ---- End of code block ----