
# Description
# To RUN: This is an independent script.
# saving_data: YES
# data_name: 
# This is to calculate the Upset Plot and volcanoplot
# Reading Materials:
# venn diagram: https://www.datanovia.com/en/blog/beautiful-ggplot-venn-diagram-with-r/
# upSet Plot: https://krassowski.github.io/complex-upset/articles/Examples_R.html


# 1. Loading Libraries ----------------------------------------------------
pacman::p_load(xlsx, tidyverse, ggplot2, ggvenn, UpSetR, ggrepel, pheatmap, glue, grDevices, svglite,tidylog)
# library(extrafont)
# loadfonts(device = "win")

# 2. Loading the data -----------------------------------------------------

source("scripts/00_path_variables.R")
factorial_anova_result <- read.csv(save_factorial_anova_result) # for venn diagram
factorial_anova_posthoc_result <- read.csv(save_factorial_anova_posthoc_result) # for upset plot
load(fold_change_result)
load(gene_data_02)
load(sample_details_01)
load(metabo_data)


# 6. Volcano Plot ---------------

source(getFoldChange)
source(getVolcanoPlot)


factorial_anova_result <- factorial_anova_result %>% 
  filter(!(Slaughter_Condition<0.05 & Feeding_Regime.Slaughter_Condition<0.05)) #removing slaughter tht shared with 

factorial_anova_result <- gene_data_02[, 1:2] %>% 
  right_join(factorial_anova_result, by = "Accession") %>% 
  left_join(fold_change_result, by = "Accession")


mrlimit = 1.2
pValue_Limit= 0.05

Factor = "Slaughter_Condition"
meanFactorLevel.1 = "Stress"
meanFactorLevel.2 = "NoStress"

columns = c("Accession", "gene_name", Factor, meanFactorLevel.1, meanFactorLevel.2)
FOldChangeData = c()
FOldChangeData = getFoldChange(factorial_anova_result, columns)

# all_FOldChangeData <- FOldChangeData %>% 
#   full_join(factorial_anova_result[, c("Accession", "Feeding_Regime.Slaughter_Condition")], by = "Accession") #%>% 
#   # filter(!Slaughter_Condition<0.05 & Feeding_Regime.Slaughter_Condition<0.05)
  

sc_volcano = getVolcanoPlot(Factor, meanFactorLevel.1, meanFactorLevel.2, FOldChangeData, mrlimit, pValue_Limit, pAdjust = TRUE, legend_title = "Slaughter Condition")
sc_volcano

ggsave(sc_volcano, file = "plots/07_05_01_volcanoplot_slaughter_condition_removing_slaughter_interaction.png", dpi = 600)

# 7. upStream and DownStream Regulator ----------------------
# Up and Down Regulated Genes
upregulated_genes <- FOldChangeData %>% 
  filter(Slaughter_Condition<0.05, log2FC>log2(1.2))

downregulated_gene <- FOldChangeData %>% 
  filter(Slaughter_Condition<0.05,log2FC< -log2(1.2))

write.xlsx(upregulated_genes, file = "data/07_05_01_up_and_down_regulated_genes_removing_interaction_terms.xlsx", sheetName = "upregulated_genes", append = TRUE)
write.xlsx(downregulated_gene, file = "data/07_05_01_up_and_down_regulated_genes_removing_interaction_terms.xlsx", sheetName = "downregulated_gene", append = TRUE)
