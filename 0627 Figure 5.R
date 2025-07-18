library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readxl)
library(purrr)
library(stringr)
library(tools)
library(patchwork)
library(tidyr)

# read CC files
file_list <- list.files(
  "D:/personal/UVA/Data analysis/202412 SILAC/Fig6/CC",
  pattern = "\\.xlsx$", full.names = TRUE
)

cc_proteins <- file_list %>%
  set_names(~ tolower(file_path_sans_ext(basename(.x)))) %>%
  map(read_excel)

library(ggrepel)
plot_cc_bubble_WT <- function(df, df_name,treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)

  light_col  <- paste0("light_", treatment, "_vs_con")
  heavy_col  <- paste0("heavy_", treatment, "_vs_con")
  total_col  <- paste0("total_", treatment, "_vs_con")
  hl_col     <- paste0("HL_", treatment, "_vs_con")
  
  process_subset <- function(cc_name, df) {
    subset_prots <- cc_proteins[[cc_name]]$Entry
    
    df %>%
      filter(PG.ProteinAccessions %in% subset_prots) %>%
      summarise(
        n_proteins = n(),
        light_mean   = mean(.data[[light_col]], na.rm = TRUE),
        light_se     = sd(.data[[light_col]], na.rm = TRUE) / sqrt(n()),
        light_median = median(.data[[light_col]], na.rm = TRUE),
        
        heavy_mean   = mean(.data[[heavy_col]], na.rm = TRUE),
        heavy_se     = sd(.data[[heavy_col]], na.rm = TRUE) / sqrt(n()),
        heavy_median = median(.data[[heavy_col]], na.rm = TRUE),
        
        total_mean   = mean(.data[[total_col]], na.rm = TRUE),
        total_se     = sd(.data[[total_col]], na.rm = TRUE) / sqrt(n()),
        total_median = median(.data[[total_col]], na.rm = TRUE),
        
        hl_mean      = mean(.data[[hl_col]], na.rm = TRUE),
        hl_se        = sd(.data[[hl_col]], na.rm = TRUE) / sqrt(n()),
        hl_median    = median(.data[[hl_col]], na.rm = TRUE)
      ) %>%
      mutate(CC = str_to_title(cc_name))
  }
  
  all_summary <- map_dfr(names(cc_proteins), ~process_subset(.x, df))
  
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
  Mito_df$HL_bin <- cut(
    Mito_df[[hl_col]],
    breaks = c(-Inf, -5, -2.5, 0, Inf),
    labels = c("H/L < -5", "-5 <= H/L < -2.5", "-2.5 <= H/L < 0", "H/L >= 0"),
    right = FALSE  # or TRUE depending on your boundary preference
  )
  p <-ggplot(all_summary, aes(x = heavy_median, y = light_median)) +
    geom_errorbarh(aes(xmin = heavy_median - heavy_se, xmax = heavy_median + heavy_se),
                   height = 0.005, color = "gray40") +
    geom_errorbar(aes(ymin = light_median - light_se, ymax = light_median + light_se),
                  width = 0.005, color = "gray40") +
    geom_point(aes(size = hl_median, color = total_median), alpha = 0.8) +
    
    geom_text(aes(label = n_proteins), vjust = 0.5, hjust = 0.5, size = 3) +
    ggrepel::geom_text_repel(
      data = all_summary %>%
        filter(!is.na(light_median), !is.na(heavy_median)),
      aes(label = CC, x = heavy_median, y = light_median),
      size = 3,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "gray60"
    )+
    scale_x_continuous(limits = c(-0.75, 0)) +
    scale_y_continuous(limits = c(0.2, 0.75)) +
    scale_size(range = c(4, 12), limits = c(-1.3, -0.5),name = "H/L ratio") +
    scale_color_gradient(low = "#6baed6", high = "#cb181d",
                         name = bquote("Total log"[2]*"(FC)"),limits = c(0,0.4), oob = scales::squish) +
    ggtitle(title_expr) +
    labs(
      x = expression("Heavy log"[2]*"(FC)"),
      y = expression("Light log"[2]*"(FC)")
    ) +
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

  ggsave("bubble_plot.pdf", plot = p, width = 6, height = 5)

  return(p)
}
plot_cc_bubble_WT(df3_cv,df_name = "df3_cv", treatment = "thg")

plot_cc_bubble_H730 <- function(df, df_name,treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)

  light_col  <- paste0("light_", treatment, "_vs_con")
  heavy_col  <- paste0("heavy_", treatment, "_vs_con")
  total_col  <- paste0("total_", treatment, "_vs_con")
  hl_col     <- paste0("HL_", treatment, "_vs_con")
  
  process_subset <- function(cc_name, df) {
    subset_prots <- cc_proteins[[cc_name]]$Entry
    
    df %>%
      filter(PG.ProteinAccessions %in% subset_prots) %>%
      summarise(
        n_proteins = n(),
        light_mean   = mean(.data[[light_col]], na.rm = TRUE),
        light_se     = sd(.data[[light_col]], na.rm = TRUE) / sqrt(n()),
        light_median = median(.data[[light_col]], na.rm = TRUE),
        
        heavy_mean   = mean(.data[[heavy_col]], na.rm = TRUE),
        heavy_se     = sd(.data[[heavy_col]], na.rm = TRUE) / sqrt(n()),
        heavy_median = median(.data[[heavy_col]], na.rm = TRUE),
        
        total_mean   = mean(.data[[total_col]], na.rm = TRUE),
        total_se     = sd(.data[[total_col]], na.rm = TRUE) / sqrt(n()),
        total_median = median(.data[[total_col]], na.rm = TRUE),
        
        hl_mean      = mean(.data[[hl_col]], na.rm = TRUE),
        hl_se        = sd(.data[[hl_col]], na.rm = TRUE) / sqrt(n()),
        hl_median    = median(.data[[hl_col]], na.rm = TRUE)
      ) %>%
      mutate(CC = str_to_title(cc_name))
  }

  all_summary <- map_dfr(names(cc_proteins), ~process_subset(.x, df))
  
  df_name_map <- list(
    "df3_cv" = "WT",
    "df4_cv" = expression("PERK"^"-/-")
  )
  
  treat_label <- ifelse(treatment == "tun", "Tm", "Th")
  
  title_expr <- if (df_name == "df4_cv") {
    bquote(PERK^-"/-" ~ .(treat_label) * "/DMSO")
  } else {
    bquote(.(df_name) ~ .(treat_label) * "/DMSO")
  }

  ggplot(all_summary, aes(x = heavy_median, y = light_median)) +
    geom_errorbarh(aes(xmin = heavy_median - heavy_se, xmax = heavy_median + heavy_se),
                   height = 0.005, color = "gray40") +
    geom_errorbar(aes(ymin = light_median - light_se, ymax = light_median + light_se),
                  width = 0.005, color = "gray40") +
    geom_point(aes(size = hl_median, color = total_median), alpha = 0.8) +
    geom_text(aes(label = n_proteins), vjust = 0.5, hjust = 0.5, size = 3) +
    ggrepel::geom_text_repel(
      data = all_summary %>%
        filter(!is.na(light_median), !is.na(heavy_median)),
      aes(label = CC, x = heavy_median, y = light_median),
      size = 3,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = Inf,
      segment.color = "gray60"
    )+
    scale_x_continuous(limits = c(-0.75, -0)) +
    scale_y_continuous(limits = c(0.2, 0.75)) +
    scale_size(range = c(4, 12), limits = c(-1.3, -0.5),name = "H/L ratio") +
    scale_color_gradient(low = "#6baed6", high = "#cb181d",
                         name = bquote("Total log"[2]*"(FC)"),limits = c(0,0.4), oob = scales::squish) +
    ggtitle(title_expr) +
    labs(
      x = expression("Heavy log"[2]*"(FC)"),
      y = expression("Light log"[2]*"(FC)")
    ) +
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
}
plot_cc_bubble_H730(df4_cv,df_name = "df4_cv", treatment = "thg")


#ATAD3A single protein Q9NVI7
plot_protein_expression_facet <- function(df_list, df_names, protein_id) {
  combined_df <- map2_dfr(df_list, df_names, function(df, name) {
    df %>%
      filter(grepl(protein_id, PG.ProteinAccessions)) %>%
      mutate(source = name)
  })
  
  gene_name <- combined_df$PG.Genes[1]
  
  protein_long <- combined_df %>%
    select(source,
           con_total_mean, tun_total_mean, thg_total_mean,
           con_total_SE, tun_total_SE, thg_total_SE,
           total_tun_vs_con_p, total_thg_vs_con_p) %>%
    pivot_longer(
      cols = -c(source, total_tun_vs_con_p, total_thg_vs_con_p),
      names_to = c("condition", "stat"),
      names_pattern = "(con|tun|thg)_total_(mean|SE)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(
      condition = factor(condition, levels = c("con", "tun", "thg"), labels = c("DMSO", "Tm", "Th")),
      channel = "Total"
    )
  
  protein_long <- protein_long %>%
    mutate(
      p_label = case_when(
        condition == "Tm" & total_tun_vs_con_p < 0.001 ~ "***",
        condition == "Tm" & total_tun_vs_con_p < 0.01 ~ "**",
        condition == "Tm" & total_tun_vs_con_p < 0.05 ~ "*",
        condition == "Th" & total_thg_vs_con_p < 0.001 ~ "***",
        condition == "Th" & total_thg_vs_con_p < 0.01 ~ "**",
        condition == "Th" & total_thg_vs_con_p < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  p <- ggplot(protein_long, aes(x = condition, y = mean, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE),
                  position = position_dodge(width = 0.8), width = 0.1) +
    geom_text(aes(label = p_label, y = mean + SE + max(mean) * 0.05),
              position = position_dodge(width = 0.8), size = 4, vjust = 0) +
    ylim(0, max(protein_long$mean + protein_long$SE, na.rm = TRUE) * 1.2) +
    facet_wrap(~source) +
    labs(
      title = gene_name,
      y = "Intensity"
    ) +
    scale_fill_manual(values = treat_colors) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      panel.spacing = unit(1, "cm")
    )
  
  return(p)
}
p_facet <- plot_protein_expression_facet(
  df_list = list(df3_cv, df4_cv),
  df_names = c("WT", "PERK-/-"),
  protein_id = "Q9NVI7"
)
print(p_facet)

#WT AND PERK KO CC comparison
df_cc_combined <- map_dfr(names(cc_proteins), function(cc_name) {
  prot_list <- cc_proteins[[cc_name]]$Entry
  
  df3_cv %>%
    filter(PG.ProteinAccessions %in% prot_list) %>%
    mutate(CC = recode(str_to_title(cc_name),
                       "Pm" = "PM",
                       "Er" = "ER",
                       "Mito" = "Mito"))
})
df_cc_combine_730 <- map_dfr(names(cc_proteins), function(cc_name) {
  prot_list <- cc_proteins[[cc_name]]$Entry
  
  df4_cv %>%
    filter(PG.ProteinAccessions %in% prot_list) %>%
    mutate(CC = recode(str_to_title(cc_name),
                       "Pm" = "PM",
                       "Er" = "ER",
                       "Mito" = "Mito"))
})
df_cc_combined <- df_cc_combined %>% mutate(Group = "WT")
df_cc_combine_730 <- df_cc_combine_730 %>% mutate(Group = "730")

df_cc_all <- bind_rows(df_cc_combined, df_cc_combine_730)

df_long <- df_cc_all %>%
  pivot_longer(
    cols = c(light_thg_vs_con, heavy_thg_vs_con, total_thg_vs_con, HL_thg_vs_con,
             light_tun_vs_con, heavy_tun_vs_con, total_tun_vs_con, HL_tun_vs_con),
    names_to = "Type",
    values_to = "log2FC"
  ) %>%
  mutate(

    Treatment = if_else(str_detect(Type, "_thg_"), "Thapsigargin", "Tunicamycin"),
    Type = case_when(
      str_detect(Type, "light_") ~ "Light",
      str_detect(Type, "heavy_") ~ "Heavy",
      str_detect(Type, "total_") ~ "Total",
      str_detect(Type, "HL_")    ~ "H/L"
    ),
    Type = factor(Type, levels = c("Light", "Heavy", "Total", "H/L")),
    Treatment = factor(Treatment, levels = c("Tunicamycin", "Thapsigargin")),
    CC = factor(CC, levels = c("Chromatin","Cytoplasm","Cytosol","Peroxisome","Proteasome","Lysosome", "Nucleus", "PM", "Golgi", "ER", "Mito", "Cytoribosome"))
  )

library(ggpubr)

#t test paired
plot_box_by_type_compare_pair <- function(df_long, type_label, treatment_label) {

  df_subset <- df_long %>%
    filter(Type == type_label, Treatment == treatment_label)

  df_paired <- df_subset %>%
    filter(!is.na(log2FC)) %>%
    select(PG.ProteinAccessions, CC, Group, log2FC) %>%
    pivot_wider(names_from = Group, values_from = log2FC) %>%
    filter(!is.na(WT), !is.na(`730`))

  signif_data <- df_paired %>%
    group_by(CC) %>%
    summarise(
      p_value = tryCatch(
        t.test(WT, `730`, paired = TRUE, alternative = "greater")$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        is.na(p_value) ~ "ns",
        p_value <= 0.001 ~ "***",
        p_value <= 0.01 ~ "**",
        p_value <= 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  count_data <- df_paired %>%
    pivot_longer(cols = c(WT, `730`), names_to = "Group", values_to = "log2FC") %>%
    group_by(CC, Group) %>%
    summarise(n = n(), .groups = "drop")
  df_plot <- df_paired %>%
    pivot_longer(cols = c(WT, `730`), names_to = "Group", values_to = "log2FC")

  p <- ggplot(df_plot, aes(x = CC, y = log2FC, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.75)) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 2, color = "black",
                 position = position_dodge(0.75)) +
    geom_text(data = signif_data,
              aes(x = CC, y = 2.1, label = label),
              inherit.aes = FALSE,
              size = 4, color = "black") +
    scale_fill_manual(values = c("WT" = "#009ADE", "730" = "#Ff1f5b")) +
    scale_y_continuous(limits = c(-2.5, 2.5))+
    labs(
      title = type_label,
      x = NULL,
      y = paste0(treatment_label, " log2 FC"),
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank()
    )
  
  return(p)
}
plot_box_by_type_compare_pair(df_long, type_label = "Heavy", treatment_label = "Thapsigargin")

#Mito
Mito <- read_excel("D:/personal/UVA/Data analysis/202412 SILAC/Fig6/CC/Mito.xlsx")
Mito_WT<- df3_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)

Mito_730<- df4_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)


plot_Mito_bubble <- function(Mito_df, df_name, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  

  light_col  <- paste0("light_", treatment, "_vs_con")
  heavy_col  <- paste0("heavy_", treatment, "_vs_con")
  total_col  <- paste0("total_", treatment, "_vs_con")
  hl_col     <- paste0("HL_", treatment, "_vs_con")
  
 
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
  
  Mito_df <- Mito_df %>%
    filter(!is.na(.data[[hl_col]]))  
  
  Mito_df$HL_bin <- cut(
    Mito_df[[hl_col]],
    breaks = c(-Inf, -5, -2.5, 0, Inf),
    labels = c("H/L < -5", "-5 <= H/L < -2.5", "-2.5 <= H/L < 0", "H/L >= 0"),
    right = FALSE  # or TRUE depending on your boundary preference
  )
 
  ggplot(Mito_df, aes(x = .data[[heavy_col]], y = .data[[light_col]])) +
    geom_point(
      aes(shape = HL_bin, color = .data[[total_col]]),
      size = 2.5, alpha = 0.8
    ) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +  
    scale_color_gradient2(
      low = "darkblue",
      mid = "#f7f7f7",
      high = "darkred",
      midpoint = 0,
      name = "Total log2(FC)",
      limits = c(-2, 3),
      oob = scales::squish
    ) +
    labs(shape = "H/L ratio bin")+
    scale_x_continuous(limits = c(-2.4,3.3))+
    scale_y_continuous(limits = c(-2,6),breaks = c(-2,0,2,4))+
    labs(
      x = expression("Heavy log"[2]*"(FC)"),
      y = expression("Light log"[2]*"(FC)"),
      title = title_expr
    ) +
    
    theme_minimal() +
    theme(
      panel.grid     = element_blank(),
      panel.border   = element_rect(color = "black", fill = NA),
      axis.line      = element_line(color = "black"),
      axis.text      = element_text(size = 8),
      axis.title     = element_text(size = 8),
      legend.position= "top",
      legend.direction = "horizontal",
      plot.title     = element_text(face = "bold", size = 8, hjust = 0.5)
    )
}
p1<-plot_Mito_bubble(Mito_WT, df_name = "df3_cv", treatment = "thg")
p2<-plot_Mito_bubble(Mito_730, df_name = "df4_cv", treatment = "thg")
combined_plot <- p1 + p2+  plot_layout(ncol = 2)
combined_plot

#Mito heavy and H+L up in WT GO
Mito_h_up <- df3_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)%>%
  filter(heavy_thg_vs_con>0.58& heavy_thg_vs_con_p<0.05&total_thg_vs_con>0.58& total_thg_vs_con_p<0.05&HL_thg_vs_con>0)

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
Mito_go <- run_go_analysis(Mito_h_up)

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
      Description = factor(Description, levels = rev(Description))  # 保持 y 轴顺序
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  # ← 显式绑定 fill
    scale_fill_gradient(low = "#ffcccc", high = "red", name = "-log10(p.adjust)") +  # 红色渐变
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment",  
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(5, 30)) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),  
      axis.line = element_line(color = "black"), 
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
}
plot_go_bubble(Mito_go , top_n = 10, title_text = "Top GO Terms (BP)")

# Mito H+L and heavy DOWN in WT GO
Mito_t_down <- df3_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)%>%
  filter(total_thg_vs_con< -0.58& total_thg_vs_con_p<0.05&heavy_thg_vs_con< -0.58& heavy_thg_vs_con_p<0.05&HL_thg_vs_con<0,!(light_thg_vs_con < -0.58 & light_thg_vs_con_p < 0.05))

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
Mito_down_go <- run_go_analysis(Mito_t_down)

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
      Description = factor(Description, levels = rev(Description))  # 保持 y 轴顺序
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  # ← 显式绑定 fill
    scale_fill_gradient(low = "lightblue", high = "blue", name = "-log10(p.adjust)") +  # 红色渐变
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment",  
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(5, 100)) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"), 
      axis.line = element_line(color = "black"), 
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
}
plot_go_bubble(Mito_down_go , top_n = 10, title_text = "Top GO Terms (BP)")

# Mito total and heavy DOWN in PERK KO
Mito_t_down_730 <- df4_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)%>%
  filter(total_thg_vs_con< -0.58& total_thg_vs_con_p<0.05&heavy_thg_vs_con< -0.58& heavy_thg_vs_con_p<0.05&HL_thg_vs_con<0,!(light_thg_vs_con < -0.58 & light_thg_vs_con_p < 0.05))

#scatter plot of PERK KO synthesis downreuglated proteins
plot_mito_scatter_down <- function(Mito_WT, Mito_730, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  
  
  light_col <- paste0("light_", treatment, "_vs_con")
  light_p   <- paste0("light_", treatment, "_vs_con_p")
  heavy_col <- paste0("heavy_", treatment, "_vs_con")
  heavy_p   <- paste0("heavy_", treatment, "_vs_con_p")
  total_col <- paste0("total_", treatment, "_vs_con")
  total_p   <- paste0("total_", treatment, "_vs_con_p")
  hl_col <- paste0("HL_", treatment, "_vs_con")
  
  merged_df <- inner_join(
    Mito_WT %>% select(PG.ProteinAccessions, PG.Genes, all_of(c(light_col, light_p,heavy_col, heavy_p, total_col, total_p,hl_col))),
    Mito_730 %>% select(PG.ProteinAccessions, PG.Genes, all_of(c(light_col, light_p,heavy_col, heavy_p, total_col, total_p,hl_col))),
    by = "PG.ProteinAccessions",
    suffix = c("_WT", "_730")
  )
  
 
  filtered_df <- merged_df %>%
    filter(
      .data[[paste0(heavy_col, "_730")]] < -0.58,
      .data[[paste0(heavy_p, "_730")]] < 0.05,
      .data[[paste0(total_col, "_730")]] < -0.58,
      .data[[paste0(total_p, "_730")]] < 0.05,
      .data[[paste0(hl_col, "_WT")]] < 0,
      !(.data[[paste0(light_col, "_WT")]] < -0.58 & .data[[paste0(light_p, "_WT")]] < 0.05)
    )
  
  
  ggplot(filtered_df, aes_string(x = paste0(heavy_col, "_WT"), y = paste0(heavy_col, "_730"))) +
    geom_point(alpha = 0.6, size = 2, color = "#4a90e2") +
    ggrepel::geom_text_repel(aes(label = PG.Genes_WT), size = 3, max.overlaps = 50) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    labs(
      x = paste("WT", heavy_col),
      y = paste("730", heavy_col),
      title = "Mito"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}
plot_mito_scatter_down(Mito_WT, Mito_730, treatment = "thg")
#GO analysis of PERK synthesis downreuglated proteins
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
Mito_down_go_730 <- run_go_analysis(Mito_t_down_730)
View(Mito_down_go_730)
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
      Description = factor(Description, levels = rev(Description))  # 保持 y 轴顺序
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  # ← 显式绑定 fill
    scale_fill_gradient(low = "lightblue", high = "blue", name = "-log10(p.adjust)") +  # 红色渐变
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment", 
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(20, 320)) + 
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),  
      axis.line = element_line(color = "black"), 
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_blank(),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )
}
plot_go_bubble(Mito_down_go_730 , top_n = 10, title_text = "Top GO Terms (BP)")
