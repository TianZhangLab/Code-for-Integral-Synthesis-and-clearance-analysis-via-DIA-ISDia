library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readxl)
library(purrr)
library(stringr)
library(tools)
library(patchwork)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(tibble)
library(pheatmap)


# read CC files
file_list <- list.files(
  "D:/personal/UVA/Data analysis/202412 SILAC/Fig6/CC",
  pattern = "\\.xlsx$", full.names = TRUE
)

cc_proteins <- file_list %>%
  set_names(~ tolower(file_path_sans_ext(basename(.x)))) %>%
  map(read_excel)

cc_combined <- imap_dfr(cc_proteins, ~ {
  df <- .x
  comp_name <- .y
  df %>%
    mutate(
      `PG.ProteinAccessions` = as.character(`Entry`),
      cellular_compartment = comp_name
    ) %>%
    dplyr::select(`PG.ProteinAccessions`, cellular_compartment)
})

df3_cv<- read.csv("D:/personal/UVA/SILAC/rebuttal/df3_cv.csv")
df4_cv<-read.csv("D:/personal/UVA/SILAC/rebuttal/df4_cv.csv")

df3_cv_CC <- df3_cv %>%
  left_join(cc_combined, by = "PG.ProteinAccessions") %>%
  mutate(cellular_compartment = replace_na(cellular_compartment, "unknown"))%>%
  dplyr::select(PG.ProteinAccessions,PG.Genes, 
                light_thg_vs_con,heavy_thg_vs_con,total_thg_vs_con, HL_thg_vs_con,light_thg_vs_con_p,heavy_thg_vs_con_p,total_thg_vs_con_p,HL_thg_vs_con_p,
                light_tun_vs_con,heavy_tun_vs_con,total_tun_vs_con, HL_tun_vs_con,light_tun_vs_con_p,heavy_tun_vs_con_p,total_tun_vs_con_p,HL_tun_vs_con_p,
                cellular_compartment
  )
df4_cv_CC <- df4_cv %>%
  left_join(cc_combined, by = "PG.ProteinAccessions") %>%
  mutate(cellular_compartment = replace_na(cellular_compartment, "unknown"))%>%
  dplyr::select(PG.ProteinAccessions,PG.Genes, light_thg_vs_con,heavy_thg_vs_con,total_thg_vs_con, HL_thg_vs_con,light_thg_vs_con_p,
                heavy_thg_vs_con_p,
                total_thg_vs_con_p,
                HL_thg_vs_con_p,
                light_tun_vs_con,heavy_tun_vs_con,total_tun_vs_con, HL_tun_vs_con,light_tun_vs_con_p,heavy_tun_vs_con_p,total_tun_vs_con_p,HL_tun_vs_con_p,
                cellular_compartment
  )
write.csv(df3_cv_CC ,"D:/personal/UVA/SILAC/Supplementary table/ISDIA_WT.csv")
write.csv(df4_cv_CC ,"D:/personal/UVA/SILAC/Supplementary table/ISDIA_PERK KO.csv")



plot_cc_bubble_WT <- function(df, df_name, treatment = c("thg", "tun")) {
  treatment <- match.arg(treatment)
  
  light_col <- paste0("light_", treatment, "_vs_con")
  heavy_col <- paste0("heavy_", treatment, "_vs_con")
  total_col <- paste0("total_", treatment, "_vs_con")
  hl_col    <- paste0("HL_", treatment, "_vs_con")
  
  process_subset <- function(cc_name, df) {
    subset_prots <- cc_proteins[[cc_name]]$Entry
    
    df_sub <- df %>%
      filter(PG.ProteinAccessions %in% subset_prots)
    
    n_light <- sum(!is.na(df_sub[[light_col]]))
    n_heavy <- sum(!is.na(df_sub[[heavy_col]]))
    n_total <- sum(!is.na(df_sub[[total_col]]))
    n_hl    <- sum(!is.na(df_sub[[hl_col]]))
    
    tibble(
      n_proteins = nrow(df_sub),
      
      light_mean   = mean(df_sub[[light_col]], na.rm = TRUE),
      light_se     = ifelse(n_light > 1, sd(df_sub[[light_col]], na.rm = TRUE) / sqrt(n_light), NA_real_),
      light_median = median(df_sub[[light_col]], na.rm = TRUE),
      
      heavy_mean   = mean(df_sub[[heavy_col]], na.rm = TRUE),
      heavy_se     = ifelse(n_heavy > 1, sd(df_sub[[heavy_col]], na.rm = TRUE) / sqrt(n_heavy), NA_real_),
      heavy_median = median(df_sub[[heavy_col]], na.rm = TRUE),
      
      total_mean   = mean(df_sub[[total_col]], na.rm = TRUE),
      total_se     = ifelse(n_total > 1, sd(df_sub[[total_col]], na.rm = TRUE) / sqrt(n_total), NA_real_),
      total_median = median(df_sub[[total_col]], na.rm = TRUE),
      
      hl_mean      = mean(df_sub[[hl_col]], na.rm = TRUE),
      hl_se        = ifelse(n_hl > 1, sd(df_sub[[hl_col]], na.rm = TRUE) / sqrt(n_hl), NA_real_),
      hl_median    = median(df_sub[[hl_col]], na.rm = TRUE),
      
      CC = str_to_title(cc_name)
    )
  }
  
  all_summary <- map_dfr(names(cc_proteins), ~ process_subset(.x, df))
  
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
  
  p <- ggplot(all_summary, aes(x = heavy_median, y = light_median)) +
    geom_errorbar(
      aes(
        xmin = heavy_median - heavy_se,
        xmax = heavy_median + heavy_se
      ),
      orientation = "y",
      height = 0.03,
      color = "gray25"
    ) +
    geom_errorbar(
      aes(
        ymin = light_median - light_se,
        ymax = light_median + light_se
      ),
      width = 0.01,
      color = "gray25"
    ) +
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
    ) +
    scale_x_continuous(limits = c(-0.75, -0.05)) +
    scale_y_continuous(limits = c(0.2, 0.75)) +
    scale_size(
      range = c(4, 12),
      limits = c(-1.3, -0.5),
      name = "H/L ratio"
    ) +
    scale_color_gradient(
      low = "#6baed6",
      high = "#cb181d",
      name = bquote("Total log"[2] * "(FC)"),
      limits = c(0, 0.4),
      oob = squish
    ) +
    ggtitle(title_expr) +
    labs(
      x = expression("Heavy log"[2] * "(FC)"),
      y = expression("Light log"[2] * "(FC)")
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
    )
  
  ggsave("bubble_plot.pdf", plot = p, width = 6, height = 5)
  
  return(p)
}
p1<- plot_cc_bubble_WT(df3_cv,df_name = "df3_cv", treatment = "thg")
p2<- plot_cc_bubble_WT(df4_cv,df_name = "df4_cv", treatment = "thg")
p3<- plot_cc_bubble_WT(df3_cv,df_name = "df3_cv", treatment = "tun")
p4<- plot_cc_bubble_WT(df4_cv,df_name = "df4_cv", treatment = "tun")

combined_plot <- p1 + p2  
print(combined_plot)


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
    dplyr::select(PG.ProteinAccessions, CC, Group, log2FC) %>%
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
#----
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
  total_p_col <- paste0("total_", treatment, "_vs_con_p")
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
    filter(!is.na(.data[[hl_col]])) %>%
    mutate(
      HL_bin = cut(
        .data[[hl_col]],
        breaks = c(-Inf, -5, -2.5, 0, Inf),
        labels = c("H/L < -5", "-5 <= H/L < -2.5", "-2.5 <= H/L < 0", "H/L >= 0"),
        right = FALSE
      ),
      Total_group = case_when(
        .data[[total_col]] > 0.58  & .data[[total_p_col]] < 0.05 ~ "Up",
        .data[[total_col]] < -0.58 & .data[[total_p_col]] < 0.05 ~ "Down",
        TRUE ~ "Not sig"
      ),
      draw_order = case_when(
        Total_group == "Not sig" ~ 1,
        Total_group %in% c("Up", "Down") ~ 2
      )
    ) %>%
    arrange(draw_order)
  
  ggplot(Mito_df, aes(x = .data[[heavy_col]], y = .data[[light_col]])) +
    geom_point(
      aes(shape = HL_bin, color = Total_group),
      size = 2.5, alpha = 0.8
    ) +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    scale_color_manual(
      values = c(
        "Up" = "#C07161",
        "Down" = "#B9A8D1",
        "Not sig" = "grey80"
      ),
      name = "Total"
    ) +
    scale_x_continuous(limits = c(-2.4, 3.3)) +
    scale_y_continuous(limits = c(-2, 3.5), breaks = c(-2, 0, 2)) +
    labs(
      x = expression("Heavy log"[2] * "(FC)"),
      y = expression("Light log"[2] * "(FC)"),
      title = title_expr,
      shape = "H/L ratio bin"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "top",
      legend.direction = "horizontal",
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5)
    )
}
p1<-plot_Mito_bubble(Mito_WT, df_name = "df3_cv", treatment = "thg")
p2<-plot_Mito_bubble(Mito_730, df_name = "df4_cv", treatment = "thg")
combined_plot <- p1 + p2+  plot_layout(ncol = 2)
combined_plot

#Mito total up in WT GO 
Mito_t_up <- df3_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)%>%
  filter(total_thg_vs_con>0.58& total_thg_vs_con_p<0.05)

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
Mito_go <- run_go_analysis(Mito_t_up)

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
    scale_fill_gradient(low = "#e8c7bf", high = "#8F3E31", name = "-log10(p.adjust)") +  # 红色渐变
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

#Heatmap
plot_one_heatmap_per_go <- function(go_ids, Mito_go, Mito_WT, Mito_730) {
  for (go_id in go_ids) {
    
    # 1. Get genes for this GO term
    go_genes_to_label <- Mito_go %>%
      filter(ID == go_id) %>%
      pull(geneID) %>%
      strsplit("/") %>%
      unlist() %>%
      unique()
    
    if (length(go_genes_to_label) == 0) {
      message("No genes found for ", go_id)
      next
    }
    
    # 2. Extract WT and PERK data
    heatmap_df <- full_join(
      Mito_WT %>%
        select(PG.Genes, heavy_thg_vs_con, light_thg_vs_con) %>%
        rename(
          heavy_WT = heavy_thg_vs_con,
          light_WT = light_thg_vs_con
        ),
      Mito_730 %>%
        select(PG.Genes, heavy_thg_vs_con, light_thg_vs_con) %>%
        rename(
          heavy_PERK = heavy_thg_vs_con,
          light_PERK = light_thg_vs_con
        ),
      by = "PG.Genes"
    ) %>%
      filter(PG.Genes %in% go_genes_to_label) %>%
      distinct(PG.Genes, .keep_all = TRUE)
    
    if (nrow(heatmap_df) == 0) {
      message("No matching genes found in Mito_WT / Mito_730 for ", go_id)
      next
    }
    
    # 3. Make matrix
    heatmap_mat <- heatmap_df %>%
      select(PG.Genes, heavy_WT, heavy_PERK, light_WT, light_PERK) %>%
      column_to_rownames("PG.Genes") %>%
      as.matrix()
    
    # 4. Get GO description
    go_title <- Mito_go %>%
      filter(ID == go_id) %>%
      summarise(title = dplyr::first(Description)) %>%
      pull(title)
    
    # 5. Draw heatmap
    pheatmap(
      heatmap_mat,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("#b9a8d1", "white", "#C07161"))(100),
      breaks = seq(-2, 2, length.out = 101),
      border_color = NA,
      fontsize_row = 8,
      fontsize_col = 10,
      main = paste0(go_id, ": ", go_title)
    )
  }
}

go_ids <- c("GO:0006520", "GO:0140053", "GO:0045333")

plot_one_heatmap_per_go(
  go_ids = go_ids,
  Mito_go = Mito_go,
  Mito_WT = Mito_WT,
  Mito_730 = Mito_730
)

#boxplot
#----
go_ids <- c("GO:0006520", "GO:0045333", "GO:0140053")
library(dplyr)
library(tidyr)
library(ggplot2)

go_ids <- c("GO:0006520", "GO:0045333", "GO:0140053")

go_gene_map <- Mito_go %>%
  filter(ID %in% go_ids) %>%
  select(ID, Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(PG.Genes = geneID, GO_category = Description) %>%
  distinct(PG.Genes, GO_category, ID)

go_gene_map

df_go_long <- full_join(
  Mito_WT %>%
    select(PG.ProteinAccessions, PG.Genes, light_thg_vs_con, heavy_thg_vs_con) %>%
    rename(
      Light = light_thg_vs_con,
      Heavy = heavy_thg_vs_con
    ) %>%
    mutate(Group = "WT"),
  
  Mito_730 %>%
    select(PG.ProteinAccessions, PG.Genes, light_thg_vs_con, heavy_thg_vs_con) %>%
    rename(
      Light = light_thg_vs_con,
      Heavy = heavy_thg_vs_con
    ) %>%
    mutate(Group = "730"),
  
  by = c("PG.ProteinAccessions", "PG.Genes"),
  suffix = c("_WTtmp", "_730tmp")
)

df_go_long <- bind_rows(
  Mito_WT %>%
    select(PG.ProteinAccessions, PG.Genes, light_thg_vs_con, heavy_thg_vs_con) %>%
    rename(
      Light = light_thg_vs_con,
      Heavy = heavy_thg_vs_con
    ) %>%
    mutate(Group = "WT"),
  
  Mito_730 %>%
    select(PG.ProteinAccessions, PG.Genes, light_thg_vs_con, heavy_thg_vs_con) %>%
    rename(
      Light = light_thg_vs_con,
      Heavy = heavy_thg_vs_con
    ) %>%
    mutate(Group = "730")
) %>%
  pivot_longer(
    cols = c(Light, Heavy),
    names_to = "Type",
    values_to = "log2FC"
  ) %>%
  inner_join(go_gene_map, by = "PG.Genes")

df_go_long <- df_go_long %>%
  mutate(
    GO_category = factor(
      GO_category,
      levels = go_gene_map %>%
        filter(ID %in% go_ids) %>%
        distinct(ID, GO_category) %>%
        pull(GO_category)
    ),
    Group = factor(Group, levels = c("WT", "730")),
    Type = factor(Type, levels = c("Light", "Heavy"))
  )

plot_box_by_go_compare_pair <- function(df_long, type_label) {
  
  df_subset <- df_long %>%
    filter(Type == type_label)
  
  df_paired <- df_subset %>%
    filter(!is.na(log2FC)) %>%
    select(PG.Genes, GO_category, Group, log2FC) %>%
    pivot_wider(names_from = Group, values_from = log2FC) %>%
    filter(!is.na(WT), !is.na(`730`))
  
  signif_data <- df_paired %>%
    group_by(GO_category) %>%
    summarise(
      p_value = tryCatch(
        t.test(WT, `730`, paired = TRUE, alternative = "greater")$p.value,
        error = function(e) NA_real_
      ),
      y_pos = max(c(WT, `730`), na.rm = TRUE) - 0.05,
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
  
  ttest_results <- df_paired %>%
    group_by(GO_category) %>%
    summarise(
      WT_mean = mean(WT, na.rm = TRUE),
      PERK_mean = mean(`730`, na.rm = TRUE),
      p_value = tryCatch(
        t.test(WT, `730`, paired = TRUE, alternative = "greater")$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  print(ttest_results)
  
  df_plot <- df_paired %>%
    pivot_longer(cols = c(WT, `730`), names_to = "Group", values_to = "log2FC")
  
  p <- ggplot(df_plot, aes(x = GO_category, y = log2FC, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.75)) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 18,
      size = 2,
      color = "black",
      position = position_dodge(0.75)
    ) +
    geom_text(
      data = signif_data,
      aes(x = GO_category, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 4,
      color = "black"
    ) +
    scale_y_continuous(limits = c(-1, 2)) +
    scale_fill_manual(values = c("WT" = "#009ADE", "730" = "#Ff1f5b")) +
    labs(
      title = type_label,
      x = NULL,
      y = "thg log2 FC",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 9, angle = 25, hjust = 1),
      axis.text.y = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank()
    )
  
  return(p)
}
p_light_go <- plot_box_by_go_compare_pair(df_go_long, "Light")
p_heavy_go <- plot_box_by_go_compare_pair(df_go_long, "Heavy")


p_light_go
p_heavy_go
p<- p_heavy_go/p_light_go
p
#----


# Mito total DOWN GO in PERK KO
Mito_t_down <- df4_cv%>%
  filter(PG.ProteinAccessions%in% Mito $Entry)%>%
  filter(total_thg_vs_con< -0.58& total_thg_vs_con_p<0.05)

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
      Description = factor(Description, levels = rev(Description))  #
    )
  
  ggplot(go_top, aes(x = FoldEnrichment, y = Description, size = Count)) +
    geom_point(aes(fill = log_p), shape = 21, color = "black", alpha = 0.8) +  
    scale_fill_gradient(low = "#e3dbef", high = "#7E68A8", name = "-log10(p.adjust)") + 
    scale_size_continuous(name = "Gene Count") +
    labs(
      title = title_text,
      x = "Fold Enrichment",  
      y = "GO Term"
    ) +
    scale_x_continuous(limits = c(0, 220),breaks = c(0,5,10)) + 
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


#heatmap 
#----
# 1. Set multiple GO IDs
go_ids <- c("GO:0006626")

# 2. Get genes from all selected GO terms
go_genes_to_label <- Mito_down_go %>%
  filter(ID %in% go_ids) %>%
  pull(geneID) %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

# 3. Extract WT and PERK data
heatmap_df <- full_join(
  Mito_WT %>%
    select(PG.Genes, heavy_thg_vs_con, light_thg_vs_con) %>%
    rename(
      heavy_WT = heavy_thg_vs_con,
      light_WT = light_thg_vs_con
    ),
  Mito_730 %>%
    select(PG.Genes, heavy_thg_vs_con, light_thg_vs_con) %>%
    rename(
      heavy_PERK = heavy_thg_vs_con,
      light_PERK = light_thg_vs_con
    ),
  by = "PG.Genes"
) %>%
  filter(PG.Genes %in% go_genes_to_label) %>%
  distinct(PG.Genes, .keep_all = TRUE)

# 4. Make matrix
heatmap_mat <- heatmap_df %>%
  select(PG.Genes, heavy_WT, heavy_PERK, light_WT, light_PERK) %>%
  column_to_rownames("PG.Genes") %>%
  as.matrix()

# 5. Make title
go_title <- Mito_down_go %>%
  filter(ID %in% go_ids) %>%
  distinct(ID, Description) %>%
  mutate(label = paste0(Description, " (", ID, ")")) %>%
  pull(label) %>%
  paste(collapse = "; ")

# 6. Draw heatmap
pheatmap(
  heatmap_mat,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#b9a8d1", "white", "#C07161"))(100),
  breaks = seq(-2, 2, length.out = 101),
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 10,
  main = go_title
)
#----