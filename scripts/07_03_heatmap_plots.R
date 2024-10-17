
# Does LogFC make sense because, it is considering the variation due to different diets,
# so volcano plots might not be the most accurate. 


# 1. Loading Libraries ----------------------------------------------------
pacman::p_load(xlsx, tidyverse, ggplot2, ggvenn, UpSetR, ggrepel, pheatmap, glue, grDevices, svglite,tidylog)

# 2. Loading the data -----------------------------------------------------

source("scripts/00_path_variables.R")
factorial_anova_result <- read.csv(save_factorial_anova_result) # for venn diagram
factorial_anova_posthoc_result <- read.csv(save_factorial_anova_posthoc_result) # for upset plot
load(fold_change_result)
load(sample_details_01)
load(metabo_data)
load(gene_data_02)

source(getFoldChange)
source(getVolcanoPlot)

# 3. Heatmap for Slaughter Conditions -------------------------------

## 3.1. Heatmap with the Residual data from model (after removing variance of feeding regime from the data) -----------------------------------
# selecting significant genes
sig_Accession_1 <- factorial_anova_result %>% 
  filter(Slaughter_Condition<0.05) %>% 
  pull(Accession)

# extracting heatmap data on the residual
heatmap_data <- read.csv(residual_model_feeding_regime_metabo) %>% 
  select(sample, all_of(sig_Accession_1))

# ordering the "sample" columns based on the "Slaughter_Condition"
sample_details_01 <- sample_details_01 %>% 
  arrange(Slaughter_Condition)

# matching the order of the rows according to ordered sample_details_01
heatmap_data <- heatmap_data %>% 
  column_to_rownames(var = "sample") %>%
  .[match(sample_details_01$sample, rownames(.)),] %>% 
  data.matrix()

# assigning annotation column
ann_col <-  factorial_anova_result %>% 
  column_to_rownames(var = "Accession") %>% 
  select(gene_name) %>% #Protein_name for plotting proteins
  rename(`gene name` = gene_name)

# assigning annotation rows
ann_row <- sample_details_01 %>% 
  column_to_rownames(var = "sample") %>% 
  select(Slaughter_Condition, Feeding_Regime) %>%  # Feeding_Regime
  rename(`Slaughter Condition` = Slaughter_Condition,
         `Feeding Regime` = Feeding_Regime)

# assigning color to different levels of Slaughter_Condition
ann_colors <-  list(
  `Slaughter Condition` = c(
    NoStress ="#FF7F0E",
    Stress = "#1F77B4"
  )
)

# plotting heatmap
heatmap_plot <- pheatmap::pheatmap(heatmap_data, 
                                   color = colorRampPalette(c("blue", "white", "red"))(100),
                                   cluster_rows = FALSE,
                                   scale = "column",
                                   annotation_row = ann_row,
                                   # annotation_col = ann_col,
                                   labels_col = ann_col$`gene name`,
                                   annotation_colors = ann_colors,
                                   show_rownames = FALSE,
                                   # annotation_legend = FALSE,
                                   angle_col = 90,
                                   gaps_row = 24,
)

png(save_heatmap_significant_genes_Slaughter_Stress_removing_affect_of_factor, width = 20, height = 10, units = 'in', res = 600) 
heatmap_plot
dev.off()


## 3.2. Heatmap for significant slaughter conditions and metabo_data ---------------

# selecting significant genes
sig_Accession_2 <- factorial_anova_result %>% 
  filter(Slaughter_Condition<0.05) %>% 
  pull(Accession)

# extracting heatmap data on the residual
heatmap_data <- metabo_data %>% 
  select(sample, all_of(sig_Accession_2))

# ordering the "sample" columns based on the "Slaughter_Condition"
sample_details_01 <- sample_details_01 %>% 
  arrange(Slaughter_Condition)

# matching the order of the rows according to ordered sample_details_01
heatmap_data <- heatmap_data %>% 
  column_to_rownames(var = "sample") %>%
  .[match(sample_details_01$sample, rownames(.)),] %>% 
  data.matrix()

# assigning annotation column
ann_col <-  factorial_anova_result %>% 
  column_to_rownames(var = "Accession") %>% 
  select(gene_name) %>% #Protein_name for plotting proteins
  rename(`gene name` = gene_name)

# assigning annotation rows
ann_row <- sample_details_01 %>% 
  column_to_rownames(var = "sample") %>% 
  select(Slaughter_Condition, Feeding_Regime) %>%  # Feeding_Regime
  rename(`Slaughter Condition` = Slaughter_Condition,
         `Feeding Regime` = Feeding_Regime)

# assigning color to different levels of Slaughter_Condition
ann_colors <-  list(
  `Slaughter Condition` = c(
    NoStress ="#FF7F0E",
    Stress = "#1F77B4"
    
  )
)
# plotting heatmap
heatmap_plot <- pheatmap::pheatmap(heatmap_data, 
                                   color = colorRampPalette(c("blue", "white", "red"))(100),
                                   cluster_rows = FALSE,
                                   scale = "column",
                                   annotation_row = ann_row,
                                   # annotation_col = ann_col,
                                   labels_col = ann_col$`gene name`,
                                   annotation_colors = ann_colors,
                                   show_rownames = FALSE,
                                   annotation_legend = TRUE,
                                   angle_col = 90,
                                   gaps_row = 24,
)

png(save_heatmap_significant_genes_Slaughter_Stress_with_metabo_data, width = 20, height = 10, units = 'in', res = 600) 
heatmap_plot
dev.off()


## 3.3. Heatmap for metabo data ad, up and down regulated slaughter conditions -------------------

IPA_Data_dest <- read.csv(IPA_Data_dest)

sig_Accession_3 <- IPA_Data_dest %>% 
  filter(p_Slaughter_Condition<0.05, 
         (log2FC_Slaughter_Condition>log2(1.2) | log2FC_Slaughter_Condition< -log2(1.2))) %>% 
  pull(Accession)

# extracting heatmap data on the residual
heatmap_data <- metabo_data %>% 
  select(sample, all_of(sig_Accession_3))

# ordering the "sample" columns based on the "Slaughter_Condition"
sample_details_01 <- sample_details_01 %>% 
  arrange(Slaughter_Condition)

# matching the order of the rows according to ordered sample_details_01
heatmap_data <- heatmap_data %>% 
  column_to_rownames(var = "sample") %>%
  .[match(sample_details_01$sample, rownames(.)),] %>% 
  data.matrix()

# assigning annotation column
ann_col <-  factorial_anova_result %>% 
  column_to_rownames(var = "Accession") %>% 
  select(gene_name) %>% #Protein_name for plotting proteins
  rename(`gene name` = gene_name)

# assigning annotation rows
ann_row <- sample_details_01 %>% 
  column_to_rownames(var = "sample") %>% 
  select(Slaughter_Condition, Feeding_Regime) %>%  # Feeding_Regime
  rename(`Slaughter Condition` = Slaughter_Condition,
         `Feeding Regime` = Feeding_Regime)

# assigning color to different levels of Slaughter_Condition
ann_colors <-  list(
  `Slaughter Condition` = c(
    NoStress ="#FF7F0E",
    Stress = "#1F77B4"
    
  )
)
# plotting heatmap
heatmap_plot <- pheatmap::pheatmap(heatmap_data, 
                                   color = colorRampPalette(c("blue", "white", "red"))(100),
                                   cluster_rows = FALSE,
                                   scale = "column",
                                   annotation_row = ann_row,
                                   # annotation_col = ann_col,
                                   labels_col = ann_col$`gene name`,
                                   annotation_colors = ann_colors,
                                   show_rownames = FALSE,
                                   annotation_legend = TRUE,
                                   angle_col = 90,
                                   gaps_row = 24,
)


png(save_heatmap_up_and_down_regulated_genes_Slaughter_Stress_with_metabo_data, width = 20, height = 10, units = 'in', res = 600) 
heatmap_plot
dev.off()



# 4. Heatmap for Feeding Regime ----------------------------------------

# selecting significant genes
sig_Accession_2 <- factorial_anova_result %>% 
  filter(Feeding_Regime<0.05) %>% 
  pull(Accession)

# extracting heatmap data on the residual
heatmap_data <- metabo_data %>% 
  select(sample, all_of(sig_Accession_2))

# ordering the "sample" columns based on the "Slaughter_Condition"
sample_details_01 <- sample_details_01 %>% 
  arrange(Slaughter_Condition)

# matching the order of the rows according to ordered sample_details_01
heatmap_data <- heatmap_data %>% 
  column_to_rownames(var = "sample") %>%
  .[match(sample_details_01$sample, rownames(.)),] %>% 
  data.matrix()

# assigning annotation column
ann_col <-  factorial_anova_result %>% 
  column_to_rownames(var = "Accession") %>% 
  select(gene_name) %>% #Protein_name for plotting proteins
  rename(`gene name` = gene_name)

# assigning annotation rows
ann_row <- sample_details_01 %>% 
  column_to_rownames(var = "sample") %>% 
  select(Slaughter_Condition, Feeding_Regime) %>%  # Feeding_Regime
  rename(`Slaughter Condition` = Slaughter_Condition,
         `Feeding Regime` = Feeding_Regime)

# assigning color to different levels of Slaughter_Condition
ann_colors <-  list(
  `Slaughter Condition` = c(
    NoStress ="#FF7F0E",
    Stress = "#1F77B4"
    
  )
)
# plotting heatmap
heatmap_plot <- pheatmap::pheatmap(heatmap_data, 
                                   color = colorRampPalette(c("blue", "white", "red"))(100),
                                   cluster_rows = FALSE,
                                   scale = "column",
                                   annotation_row = ann_row,
                                   # annotation_col = ann_col,
                                   labels_col = ann_col$`gene name`,
                                   annotation_colors = ann_colors,
                                   show_rownames = FALSE,
                                   annotation_legend = TRUE,
                                   angle_col = 90,
                                   gaps_row = 24,
)

png(save_heatmap_significant_genes_Feeding_Regime_with_metabo_data, width = 20, height = 10, units = 'in', res = 600) 
heatmap_plot
dev.off()