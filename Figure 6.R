library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

complex_results_th <- read_excel("complex_results.xlsx",
                                 sheet = "Complex Th intensity results")

##Figure 6A
# Extract Th data (Heavy and Light) with q-values
# Extract Th data (Heavy and Light) with q-values
th_df <- complex_results_th[, c("complex_name", "ThHeavy_WT_effect", "ThHeavy_WT_se", "ThHeavy_WT_qval",
                                "ThLight_WT_effect", "ThLight_WT_se", "ThLight_WT_qval")]
colnames(th_df) <- c("complex_name", "ThHeavy_WT", "ThHeavy_WT_se", "ThHeavy_qval",
                     "ThLight_WT", "ThLight_WT_se", "ThLight_qval")

# Define significance categories
th_df$Significance <- "None"
th_df$Significance[th_df$ThHeavy_qval < 0.05 & th_df$ThLight_qval < 0.05] <- "Both"
th_df$Significance[th_df$ThHeavy_qval < 0.05 & th_df$ThLight_qval >= 0.05] <- "Heavy only"
th_df$Significance[th_df$ThHeavy_qval >= 0.05 & th_df$ThLight_qval < 0.05] <- "Light only"

# Function to select top 10 by ThHeavy_WT in each category
get_top10 <- function(df, category) {
  df %>%
    filter(Significance == category) %>%
    arrange(-ThHeavy_WT) %>%
    slice(1:10)
}

# Combine top 10 from each category
top10_both <- get_top10(th_df, "Both")
top10_heavy <- get_top10(th_df, "Heavy only")
top10_light <- get_top10(th_df, "Light only")
top10_none <- get_top10(th_df, "None")

# Combine all top 10s
top10_df <- bind_rows(top10_both, top10_heavy, top10_light, top10_none)

# Subset to all "Light only" points
light_only_df <- th_df %>%
  filter(Significance == "Light only")

# Combine "Light only" and "top 10" labels
label_df <- bind_rows(light_only_df, top10_df) %>%
  distinct(complex_name, .keep_all = TRUE)  # Remove duplicates if any

top_labels <- th_df %>%
  filter(complex_name == "Spliceosome-Sm") 

ggplot(th_df, aes(x = ThHeavy_WT, y = ThLight_WT, color = Significance)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ThLight_WT - ThLight_WT_se, ymax = ThLight_WT + ThLight_WT_se), width = 0.02) +
  geom_errorbarh(aes(xmin = ThHeavy_WT - ThHeavy_WT_se, xmax = ThHeavy_WT + ThHeavy_WT_se), height = 0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +   # Make x=0 darker
  geom_hline(yintercept = 0, color = "black", size = 0.5) +   # Make y=0 darker
  geom_text(data = label_df, aes(label = complex_name), vjust = -1, hjust = 0.5, size = 3, color = "black") +
  #geom_text(data = top_labels, aes(label = complex_name), vjust = -1, size = 3, color = "black") +
  scale_color_manual(values = c("Both" = "orange", "Heavy only" = "red", "Light only" = "purple", "None" = "gray")) +
  labs(
    title = "Th Condition: Heavy vs Light with Labels for 'Light only' & Top 10 per Category",
    x = "Th Heavy WT Effect",
    y = "Th Light WT Effect",
    color = "Significance"
  ) +
  theme_minimal()

# Figure 6B

df <- read_excel("C:/Users/cgm8ck/Downloads/complex_results.xlsx")

# Set significance threshold
qval_threshold <- 0.05

##heavy
# Create a significance category
df <- df %>%
  mutate(Significance = case_when(
    ThHeavy_WT_qval < qval_threshold & ThHeavy_KO_qval < qval_threshold ~ "Significant in Both",
    ThHeavy_WT_qval < qval_threshold & ThHeavy_KO_qval >= qval_threshold ~ "Significant in WT only",
    ThHeavy_WT_qval >= qval_threshold & ThHeavy_KO_qval < qval_threshold ~ "Significant in KO only",
    TRUE ~ "Not Significant"
  ))

# Find the top 2 complexes in each category by ThHeavy_WT_effect
top_labels <- df %>%
  group_by(Significance) %>%
  top_n(10, -ThHeavy_WT_effect) %>%
  ungroup()

# Define Cell-style colors
cell_colors <- c("Significant in Both" = "orange", 
                 "Significant in WT only" = "#4DBBD5",
                 "Significant in KO only" = "#00A087",
                 "Not Significant" = "#B09C85")


ggplot(df, aes(x = ThHeavy_WT_effect, y = ThHeavy_KO_effect, color = Significance, fill = Significance)) +
  geom_point(size = 3, shape = 21, stroke = 0.5) +
  geom_errorbar(aes(ymin = ThHeavy_KO_effect - ThHeavy_KO_se, ymax = ThHeavy_KO_effect + ThHeavy_KO_se), width = 0.02) +
  geom_errorbarh(aes(xmin = ThHeavy_WT_effect - ThHeavy_WT_se, xmax = ThHeavy_WT_effect + ThHeavy_WT_se), height = 0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_text(data = top_labels, aes(label = complex_name), vjust = -1, size = 3, color = "black") +
  scale_color_manual(values = cell_colors) +
  scale_fill_manual(values = cell_colors) +
  labs(
    title = "ThHeavy Effect (WT vs KO) - Top Complexes Labeled",
    x = "Th Heavy Effect in WT",
    y = "Th Heavy Effect in KO",
    color = "Significance",
    fill = "Significance"
  ) +
  theme_minimal()

# Figure 6C & Fig. S10C
library(readxl)
library(dplyr)
library(tidyr)
library(cowplot)
a012_wt_dat <- read_excel("A012_WT_data.xlsx")
h730_ko_dat <- read_excel("H730_knockout_data.xlsx")

## Filter out NA rows
a012_wt_dat <- a012_wt_dat[apply(a012_wt_dat[,-c(1,2)], 1, function(x) !all(is.na(x))),]
h730_ko_dat <- h730_ko_dat[apply(h730_ko_dat[,-c(1,2)], 1, function(x) !all(is.na(x))),]

#### Pull link from UniProt to gene symbol
a012_uniprot_to_gene <- a012_wt_dat %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes) %>%
  distinct %>%
  dplyr::rename(external_gene_name = PG.Genes)

#### Pull link from UniProt to gene symbol
h730_uniprot_to_gene <- h730_ko_dat %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes) %>%
  distinct %>%
  dplyr::rename(external_gene_name = PG.Genes)

#### Protein complexes defined in Ori et al
protein_complex_dat <- read_excel("ori_etal_protein_complexes.xlsx") %>%
  dplyr::rename(complex_name = "Complex Name",
                num_human_members = "# of members",
                human_ensembl_ids = "Member identifiers (human Ensembl gene)")

biomart_results <- readRDS("biomart_results.rds")

#### Pulling human ENSEMBL gene IDs
complex_ids <- sapply(1:nrow(protein_complex_dat), 
                      function(i) rep(protein_complex_dat$`Complex ID`[i], length(unlist(strsplit(protein_complex_dat$human_ensembl_ids[i], split = " ")))), 
                      simplify = TRUE) %>%
  unlist

human_gene_ids <- unlist(sapply(1:nrow(protein_complex_dat), 
                                function(i) strsplit(protein_complex_dat$human_ensembl_ids[i], split = " ")))

complex_to_gene <- data.frame(ensembl_gene_id = human_gene_ids,
                              complex_id = complex_ids) %>%
  left_join(protein_complex_dat %>%
              dplyr::rename(complex_id = "Complex ID") %>%
              dplyr::select(complex_id, complex_name))

complex_to_gene <- complex_to_gene %>%
  left_join(biomart_results)

## Filter out NA rows
a012_wt_dat <- a012_wt_dat[apply(a012_wt_dat[,-c(1,2)], 1, function(x) !all(is.na(x))),]
h730_ko_dat <- h730_ko_dat[apply(h730_ko_dat[,-c(1,2)], 1, function(x) !all(is.na(x))),]

#### Pull link from UniProt to gene symbol
a012_uniprot_to_gene <- a012_wt_dat %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes) %>%
  distinct %>%
  dplyr::rename(external_gene_name = PG.Genes)

#### Pull link from UniProt to gene symbol
h730_uniprot_to_gene <- h730_ko_dat %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes) %>%
  distinct %>%
  dplyr::rename(external_gene_name = PG.Genes)

a012_wt_long_dat <- a012_wt_dat %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes, 
                Th3_Light, Th2_Light, Th1_Light,
                Th3_Heavy, Th2_Heavy, Th1_Heavy,
                Th3_Total, Th2_Total, Th1_Total,
                Tm3_Light, Tm2_Light, Tm1_Light,
                Tm3_Heavy, Tm2_Heavy, Tm1_Heavy,
                Tm3_Total, Tm2_Total, Tm1_Total,
                DMSO3_Light, DMSO2_Light, DMSO1_Light,
                DMSO3_Heavy, DMSO2_Heavy, DMSO1_Heavy,
                DMSO3_Total, DMSO2_Total, DMSO1_Total) %>%
  distinct %>%
  pivot_longer(-c(PG.ProteinAccessions, PG.Genes), names_to = "sample_id", values_to = "intensity") %>%
  mutate(treatment = case_when(grepl(x = sample_id, pattern = "Th") ~ "Th",
                               grepl(x = sample_id, pattern = "Tm") ~ "Tm",
                               grepl(x = sample_id, pattern = "DMSO") ~ "DMSO"),
         experiment = case_when(grepl(x = sample_id, pattern = "Light") ~ "Light",
                                grepl(x = sample_id, pattern = "Heavy") ~ "Heavy",
                                grepl(x = sample_id, pattern = "Total") ~ "Total")) %>%
  mutate(category = paste(treatment, experiment)) %>%
  dplyr::rename(external_gene_name = PG.Genes) %>%
  mutate(experiment = factor(experiment, levels = c("Light", "Heavy", "Total"))) %>%
  left_join(complex_to_gene, relationship = "many-to-many")

a012_wt_long_dat$intensity <- as.numeric(a012_wt_long_dat$intensity)
a012_wt_dmso_mean_dat <- a012_wt_long_dat %>%
  filter(treatment == "DMSO") %>%
  group_by(PG.ProteinAccessions, experiment) %>%
  summarize(mean_dmso_intensity = mean(log2(intensity), na.rm = TRUE)) %>%
  ungroup
a012_wt_long_dat <- a012_wt_long_dat %>%
  left_join(a012_wt_dmso_mean_dat) %>%
  mutate(gene_label = ifelse(!is.na(external_gene_name), external_gene_name, paste("UniProt:", PG.ProteinAccessions, sep = "\n")))


#730
h730_ko_long_dat <- h730_ko_dat %>%
  dplyr::select(PG.ProteinAccessions, PG.Genes, 
                Th3_Light, Th2_Light, Th1_Light,
                Th3_Heavy, Th2_Heavy, Th1_Heavy,
                Th3_Total, Th2_Total, Th1_Total,
                Tm3_Light, Tm2_Light, Tm1_Light,
                Tm3_Heavy, Tm2_Heavy, Tm1_Heavy,
                Tm3_Total, Tm2_Total, Tm1_Total,
                DMSO3_Light, DMSO2_Light, DMSO1_Light,
                DMSO3_Heavy, DMSO2_Heavy, DMSO1_Heavy,
                DMSO3_Total, DMSO2_Total, DMSO1_Total) %>%
  distinct %>%
  pivot_longer(-c(PG.ProteinAccessions, PG.Genes), names_to = "sample_id", values_to = "intensity") %>%
  mutate(treatment = case_when(grepl(x = sample_id, pattern = "Th") ~ "Th",
                               grepl(x = sample_id, pattern = "Tm") ~ "Tm",
                               grepl(x = sample_id, pattern = "DMSO") ~ "DMSO"),
         experiment = case_when(grepl(x = sample_id, pattern = "Light") ~ "Light",
                                grepl(x = sample_id, pattern = "Heavy") ~ "Heavy",
                                grepl(x = sample_id, pattern = "Total") ~ "Total")) %>%
  mutate(category = paste(treatment, experiment)) %>%
  dplyr::rename(external_gene_name = PG.Genes) %>%
  mutate(experiment = factor(experiment, levels = c("Light", "Heavy", "Total"))) %>%
  left_join(complex_to_gene, relationship = "many-to-many")

h730_ko_long_dat$intensity <- as.numeric(h730_ko_long_dat$intensity)
h730_ko_dmso_mean_dat <- h730_ko_long_dat %>%
  filter(treatment == "DMSO") %>%
  group_by(PG.ProteinAccessions, experiment) %>%
  summarize(mean_dmso_intensity = mean(log2(intensity), na.rm = TRUE)) %>%
  ungroup
h730_ko_long_dat <- h730_ko_long_dat %>%
  left_join(h730_ko_dmso_mean_dat) %>%
  mutate(gene_label = ifelse(!is.na(external_gene_name), external_gene_name, paste("UniProt:", PG.ProteinAccessions, sep = "\n")))


# Define the complex of interest
complex_of_interest <- "p97/VCP-VIMP-DERL1-DERL2-HRD1-SEL1L complex"

library(dplyr)
library(ggplot2)
library(ggbeeswarm)


# Filter and compute log2 relative intensity for WT
wt_heavy <- a012_wt_long_dat %>%
  filter(category %in% c("Th Heavy", "DMSO Heavy"),
         complex_name == complex_of_interest) %>%
  mutate(log2_rel_intensity = log2(intensity) - mean_dmso_intensity,
         Group = "WT")

wt_light <- a012_wt_long_dat %>%
  filter(category %in% c("Th Light", "DMSO Light"),
         complex_name == complex_of_interest) %>%
  mutate(log2_rel_intensity = log2(intensity) - mean_dmso_intensity,
         Group = "WT")

wt_total <- a012_wt_long_dat %>%
  filter(category %in% c("Th Total", "DMSO Total"),
         complex_name == complex_of_interest) %>%
  mutate(log2_rel_intensity = log2(intensity) - mean_dmso_intensity,
         Group = "WT")

#make the plots

plight <- ggplot(wt_light, aes(x = category, y = log2_rel_intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1, alpha = 0.7) +
  facet_wrap(~PG.ProteinAccessions, scales = "free_y", ncol=7) +
  labs(
    title = paste("Heavy Comparison (WT vs KO) for", complex_of_interest),
    x = "Category",
    y = "Log2 Relative Intensity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("WT"="#ffbb6f"))

pheavy <- ggplot(wt_heavy, aes(x = category, y = log2_rel_intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1, alpha = 0.7) +
  facet_wrap(~PG.ProteinAccessions, scales = "free_y", ncol=7) +
  labs(
    title = paste("Heavy Comparison (WT vs KO) for", complex_of_interest),
    x = "Category",
    y = "Log2 Relative Intensity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("WT"="#999999"))

ptotal <- ggplot(wt_total, aes(x = category, y = log2_rel_intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1, alpha = 0.7) +
  facet_wrap(~PG.ProteinAccessions, scales = "free_y", ncol=7) +
  labs(
    title = paste("Heavy Comparison (WT vs KO) for", complex_of_interest),
    x = "Category",
    y = "Log2 Relative Intensity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("WT"="#5e4c5f"))
plot_grid(ptotal, pheavy, plight, nrow = 3)

##ko plots
ko_heavy <- h730_ko_long_dat %>%
  filter(category %in% c("Th Heavy", "DMSO Heavy"),
         complex_name == complex_of_interest) %>%
  mutate(log2_rel_intensity = log2(intensity) - mean_dmso_intensity,
         Group = "KO")

ko_light <- h730_ko_long_dat %>%
  filter(category %in% c("Th Light", "DMSO Light"),
         complex_name == complex_of_interest) %>%
  mutate(log2_rel_intensity = log2(intensity) - mean_dmso_intensity,
         Group = "KO")

ko_total <- h730_ko_long_dat %>%
  filter(category %in% c("Th Total", "DMSO Total"),
         complex_name == complex_of_interest) %>%
  mutate(log2_rel_intensity = log2(intensity) - mean_dmso_intensity,
         Group = "KO")

#make the plots

pkolight <- ggplot(ko_light, aes(x = category, y = log2_rel_intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1, alpha = 0.7) +
  facet_wrap(~PG.ProteinAccessions, scales = "free_y", ncol=7) +
  labs(
    title = paste("Heavy Comparison (WT vs KO) for", complex_of_interest),
    x = "Category",
    y = "Log2 Relative Intensity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("WT"="#ffbb6f"))

pkoheavy <- ggplot(ko_heavy, aes(x = category, y = log2_rel_intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1, alpha = 0.7) +
  facet_wrap(~PG.ProteinAccessions, scales = "free_y", ncol=7) +
  labs(
    title = paste("Heavy Comparison (WT vs KO) for", complex_of_interest),
    x = "Category",
    y = "Log2 Relative Intensity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("WT"="#999999"))

pkototal <- ggplot(ko_total, aes(x = category, y = log2_rel_intensity, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), size = 1, alpha = 0.7) +
  facet_wrap(~PG.ProteinAccessions, scales = "free_y", ncol=7) +
  labs(
    title = paste("Heavy Comparison (WT vs KO) for", complex_of_interest),
    x = "Category",
    y = "Log2 Relative Intensity"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = c("WT"="#5e4c5f"))
plot_grid(pkototal, pkoheavy, pkolight, nrow = 3)


