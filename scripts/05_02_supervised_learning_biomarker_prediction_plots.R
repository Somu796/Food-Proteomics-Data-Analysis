# This script only can be run when 06_02 is already done

# Variable Importance: https://youtu.be/O8Un6lQlnB8?si=2tDjM5n5KeTBwkEH
# vivid: https://youtu.be/noHUJeouff8?si=pqREf_b7U-o5htXF , for nice variable importance pots and their two way interactions
# machine learning in R: https://www.youtube.com/live/Liws4MShq1A?si=LlCnp0LuBFjEnT


# 1. install and import the package ---------------------------------------
pacman::p_load(tidyverse, janitor, EFS, ggplot2, ggrepel, ggpubr, ComplexHeatmap, tidylog)

source("scripts/00_path_variables.R")
source("scripts/functions/getVolcanoPlot2.R")

source(getFoldChange)
source(getVolcanoPlot)

# 2. Importing Data ----------------------------------------------

# efs_result_dest <- "data/temp_Ensemble_Feature_selection_output.csv"
important_genes_parts <- read.csv(efs_result_dest) %>% 
  filter(total >=0.45) %>% 
  column_to_rownames("gene_name")

# factorial_anova_result <- read.csv(factorial_anova_result)
load(gene_data_02)
# load(fold_change_result)
load(sample_details_01)

# factorial_anova_result_table <- gene_data_02[,c("Accession", "gene_name")] %>%
#   right_join(factorial_anova_result, by = "Accession")



# # 3. Plotting -------------------------------------------------------------
# 
# ## 3.1. volacno Plot EFS ---------------------------
# mrlimit = 1.2
# pValue_Limit= 0.05
# 
# factorial_anova_result <- gene_data_02[, 1:2] %>% 
#   right_join(factorial_anova_result, by = "Accession") %>% 
#   full_join(fold_change_result, by = "Accession")
# 
# Factor = "Slaughter_Condition"
# meanFactorLevel.1 = "NoStress"
# meanFactorLevel.2 = "Stress"
# 
# columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
# FOldChangeData = c()
# FOldChangeData = getFoldChange(factorial_anova_result, columns)
# 
# filtered_data <- FOldChangeData[FOldChangeData$Accession %in% rownames(important_genes_parts),]
# 
# volcano_plot_EFS <- getVolcanoPlotEFS(Factor, meanFactorLevel.1, meanFactorLevel.2 ,FOldChangeData, mrlimit, pValue_Limit, pAdjust = TRUE, legend_title = "Horn Status", additional_data = filtered_data)
# 
# windows(16, 10)
# volcano_plot_EFS
# 
# ggsave(volcano_plot_EFS, file = save_volcano_plot_EFS, dpi = 600)
# 
# ## 3.2. Volcano Plot EFS vs ANOVA ---------------------------------
# horn_volcano = getVolcanoPlot(Factor, meanFactorLevel.1, meanFactorLevel.2 ,FOldChangeData, mrlimit, pValue_Limit, pAdjust = TRUE, legend_title = "Horn Status") #V = 
# 
# EFS_vs_ANOVA <- ggarrange(horn_volcano,volcano_plot_EFS, nrow = 1, ncol = 2) 
# 
# windows(16,10)
# EFS_vs_ANOVA
# 
# ggsave(EFS_vs_ANOVA, file = save_EFS_vs_ANOVA, dpi = 600)

## 3.3. Heat map --------------------------

load(metabo_data)

metabo_data <- sample_details_01 %>% 
  full_join(metabo_data, by = "sample")


data_01 <- metabo_data
col_names_heat_map <- c("sample", colnames(data_01)[colnames(data_01) %in% rownames(important_genes_parts)])

heat_map_data <- data_01 %>% 
  arrange(factor(Slaughter_Condition)) 

ann_col <- heat_map_data %>%  select(Slaughter_Condition, Feeding_Regime) %>% rename("Slaughter Condition" = "Slaughter_Condition", "Feeding Regime" = "Feeding_Regime")
rownames(ann_col) <- heat_map_data$sample

df <- heat_map_data %>% 
  select(all_of(col_names_heat_map)) %>% 
  t() %>% 
  as.data.frame() %>% 
  row_to_names(1) %>% 
  data.matrix()

ann_colors = list(
  "Slaughter Condition"  = c("Stress"="#F0978D", "NoStress"= "#9cadea"))


heatmap_plot <- pheatmap(
  df, 
  cluster_cols = FALSE, 
  scale = "row",
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  # border_color = "white",
  # gaps_col = cumsum(c(sum(heat_map_data$Horns == "Horn"),sum(heat_map_data$Horns == "Horn")+1)),
  legend = TRUE,
  heatmap_legend_param = list(title = "Intensity"),
  # cutree_cols = 2
)

windows(16, 10)
heatmap_plot

png(save_heatmap_EFS_features_Slaughter_Stress, width = 16, height = 10, units = 'in', res = 600) 
heatmap_plot
dev.off()

